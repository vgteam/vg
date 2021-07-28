#include "seed_clusterer.hpp"

#include <algorithm>

//#define DEBUG_CLUSTER
namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer( MinimumDistanceIndex& dist_index) :
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

    //Helper function to get the minimum value that is not -1
    int64_t min_not_minus_one(int64_t n1, int64_t n2) {
        return static_cast<int64_t>(std::min(static_cast<uint64_t>(n1), static_cast<uint64_t>(n2)));
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

        //For each level of the snarl tree, maps snarls (index into
        //dist_index.snarl_indexes) at that level to nodes belonging to the snarl
        //This is later used to populate snarl_to_node in the tree state
        vector<hash_map<size_t, vector<pair<NetgraphNode, NodeClusters>>>> snarl_to_nodes_by_level;
        snarl_to_nodes_by_level.reserve(dist_index.tree_depth+1);



        //This stores all the tree relationships and cluster information
        //for a single level of the snarl tree as it is being processed
        //It also keeps track of the parents of the current level
        size_t seed_count = 0;
        for (auto v : all_seeds) seed_count+= v->size();
        TreeState tree_state (&all_seeds, read_distance_limit, fragment_distance_limit, seed_count);


        //Populate tree_state.node_to_seeds (mapping each node to the seeds it
        //contains) and snarl_to_nodes_by_level
        get_nodes(tree_state, snarl_to_nodes_by_level);

        //Initialize the tree state to be the bottom level
        tree_state.snarl_to_nodes = std::move(snarl_to_nodes_by_level[snarl_to_nodes_by_level.size() - 1]);

        for (int depth = snarl_to_nodes_by_level.size() - 1 ; depth >= 0 ; depth --) {
            //Go through each level of the tree, bottom up, and cluster that
            // level. Each level includes the snarl at that level, the nodes
            // belonging to those snarls, and the chains comprised of them
            //
            // tree_state knows all children of the snarls at this level

            // Bring in the direct child nodes that come in at this level in the snarl tree.
            // They only ever occur below the root.
            if (depth != 0) {
                tree_state.parent_snarl_to_nodes = std::move(snarl_to_nodes_by_level[depth - 1]);
            }

#ifdef DEBUG_CLUSTER
assert(tree_state.read_index_offsets[0] == 0);
for (size_t i = 1 ; i < tree_state.all_seeds->size() ; i++) {
    assert (tree_state.read_index_offsets[i] + tree_state.all_seeds->at(i)->size() == tree_state.read_index_offsets[i+1]);
}
#endif
            //Cluster all the snarls at this depth
            //Also records which snarls are in chains and the parents of these
            //snarls in tree_state.parent_snarl_to_node
            cluster_snarl_level(tree_state, depth);

            //And cluster all the chains, record the parents of these chains
            cluster_chain_level(tree_state, depth);

            // Swap buffer over for the next level
            tree_state.snarl_to_nodes = std::move(tree_state.parent_snarl_to_nodes);
            tree_state.chain_to_snarls.clear();
        }


#ifdef DEBUG_CLUSTER

        cerr << "Found read clusters : " << endl;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << "\t read num " << read_num << " : " << endl;
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
                        int64_t d2 = min_not_minus_one(d1, dist_index.min_distance(pos1, rev2));
                        int64_t d3 = min_not_minus_one(d2, dist_index.min_distance(rev1, rev2));
                        int64_t d4 = min_not_minus_one(d3, dist_index.min_distance(rev1, pos2));
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
                               int64_t d2 = min_not_minus_one(d1, dist_index.min_distance(pos1, rev2));
                               int64_t d3 = min_not_minus_one(d2, dist_index.min_distance(rev1, rev2));
                               int64_t d4 = min_not_minus_one(d3, dist_index.min_distance(rev1, pos2));
                               
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


    void SnarlSeedClusterer::get_nodes( TreeState& tree_state,
              vector<hash_map<size_t,vector<pair<NetgraphNode, NodeClusters>>>>&
                                                 snarl_to_nodes_by_level) const {
#ifdef DEBUG_CLUSTER
cerr << "Nested positions: " << endl << "\t";
#endif

        // Assign each seed to a node.
        hash_set<id_t> seen_nodes;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++){ 
            const vector<Seed>* seeds = tree_state.all_seeds->at(read_num);
            for (size_t i = 0; i < seeds->size(); i++) {
                pos_t pos = seeds->at(i).pos;
                id_t id = get_id(pos);
                

                //Assign the seed to a node
                if (!seeds->at(i).is_top_level_node && !seeds->at(i).is_top_level_snarl) {
                    //If this seed is not on a top-level chain or top-level simple bubble
                    //A seed can still be added here if it is on a top-level chain
                    tree_state.node_to_seeds[read_num].emplace_back(id, i);
#ifdef DEBUG_CLUSTER
                    cerr << read_num << ":" << pos << ", ";
#endif

                    //And the node to a snarl
                    if (seen_nodes.count(id) < 1) {
                        seen_nodes.insert(id);
                        size_t snarl_i = dist_index.get_primary_assignment(id);
                        size_t depth = dist_index.snarl_indexes[snarl_i].depth;
                        if (depth+1 > snarl_to_nodes_by_level.size()) {
                            snarl_to_nodes_by_level.resize(depth+1);
                        }
                        snarl_to_nodes_by_level[depth][snarl_i].emplace_back(
                                 NetgraphNode(id, NODE), NodeClusters(tree_state.all_seeds->size()));
                    } 
                } else if (seeds->at(i).is_top_level_node) {
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

                    tree_state.node_to_seeds[read_num].emplace_back(id, i);

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
            }
            std::sort(tree_state.node_to_seeds[read_num].begin(), tree_state.node_to_seeds[read_num].end());
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


    void SnarlSeedClusterer::cluster_snarl_level(TreeState& tree_state, size_t depth) const {

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

                size_t chain_assignment = dist_index.get_chain_assignment(snarl_index.parent_id);
                size_t chain_rank = dist_index.get_chain_rank(snarl_index.id_in_parent);

                tree_state.chain_to_snarls[chain_assignment].emplace(
                        chain_rank, make_pair(snarl_i, cluster_one_snarl(tree_state, snarl_i)));

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
                    size_t parent_snarl_i = dist_index.get_primary_assignment( snarl_index.parent_id);

                    tree_state.parent_snarl_to_nodes[parent_snarl_i].emplace_back(
                            NetgraphNode (snarl_i, SNARL), cluster_one_snarl(tree_state, snarl_i));

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

                    NodeClusters top_snarl_clusters = cluster_one_snarl(tree_state, snarl_i);
                }
            }
        }
    }

    void SnarlSeedClusterer::cluster_chain_level(TreeState& tree_state, size_t depth) const {
        vector<bool>  seen_components( tree_state.top_level_seed_clusters.size(), false);
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

            //Mark this component as being seen
            size_t component = dist_index.get_connected_component(dist_index.chain_indexes[chain_i].id_in_parent);
            if (tree_state.component_to_index.count(component) != 0) {
                seen_components[tree_state.component_to_index[component]] = true;
            }
            // Compute the clusters for the chain
            if (depth == 0) {
                cluster_one_chain(tree_state, chain_i, depth);
            } else {
                NodeClusters chain_clusters = cluster_one_chain(tree_state, chain_i, depth);

                // We actually have a parent

                // Find the node ID that heads the parent of that chain.
                size_t parent_id = dist_index.chain_indexes[chain_i].parent_id;
                // It must be a legitimate node ID we cover.
                assert(parent_id >= dist_index.min_node_id);
                assert(parent_id <= dist_index.max_node_id);

                // Map it to the snarl number that should be represented by it
                // (and thus also contain the chain)
                size_t parent_snarl_i = dist_index.get_primary_assignment(parent_id);

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
        if (depth == 0) {
            //If this is the top-level, go through components for which there were no nested seeds
            // and cluster the top-level seeds and snarls
            for (size_t component_num = 0 ; component_num < seen_components.size() ; component_num++) {
                if (!seen_components[component_num]) {
                    if (tree_state.simple_snarl_to_nodes_by_component[component_num].empty()) {
                        //If there are no top-level simple bubbles in this component, cluster only top level seeds
                        cluster_only_top_level_chain_seeds(tree_state, tree_state.top_level_seed_clusters[component_num]);
                    } else {
                        //Cluster both top-level bubbles and chain nodes
#ifdef DEBUG_CLUSTER
                        cerr << "Clustering top-level bubbles and nodes " << endl;
#endif
                        id_t node_id = std::get<0>(tree_state.simple_snarl_to_nodes_by_component[component_num].begin()->second.front());
                        size_t chain_i = dist_index.component_to_chain_index[dist_index.get_connected_component(node_id)-1];
                        cluster_one_chain(tree_state, chain_i, depth);
                    }
                }
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
        NodeClusters node_clusters(tree_state.all_seeds->size());

        if (tree_state.read_distance_limit > node_length) {
            //If the limit is greater than the node length, then all the
            //seeds on this node must be in the same cluster

            size_t fragment_group_id = -1;
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                auto seed_range_start = std::lower_bound(
                    tree_state.node_to_seeds[read_num].begin(), tree_state.node_to_seeds[read_num].end(),
                    std::pair<id_t, size_t>(node_id, 0));
                if (seed_range_start != tree_state.node_to_seeds[read_num].end() 
                        && seed_range_start->first == node_id) {

                    size_t group_id = seed_range_start->second;

                    for (auto iter = seed_range_start; iter != tree_state.node_to_seeds[read_num].end() 
                                                      && iter->first == node_id; ++iter) {
                        //For each seed on this node, add it to the cluster
                        //And find the shortest distance from any seed to both
                        //ends of the node

                        pos_t seed = tree_state.all_seeds->at(read_num)->at(iter->second).pos;
                        int64_t dist_left = is_rev(seed) ? node_length- get_offset(seed) : get_offset(seed) + 1;
                        int64_t dist_right = is_rev(seed) ? get_offset(seed) + 1 : node_length - get_offset(seed);

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
                            if (fragment_group_id == -1 ) {
                                fragment_group_id = seed_range_start->second + tree_state.read_index_offsets[read_num];
                            }
                            tree_state.fragment_union_find.union_groups(
                                    fragment_group_id, iter->second + tree_state.read_index_offsets[read_num]);
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
                    //Record the previous cluster
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
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
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



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_snarl(
                    TreeState& tree_state, size_t snarl_index_i) const {
        /*Get the clusters on this snarl.
         * Nodes have not yet been clustered */
        MinimumDistanceIndex::SnarlIndex& snarl_index = dist_index.snarl_indexes[snarl_index_i];
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl number " << snarl_index_i
             << " headed by node " << snarl_index.id_in_parent << " and running to " << snarl_index.end_id << endl;
        bool at_snarl = true;
        id_t at_id = snarl_index.id_in_parent;
        size_t at_rank = snarl_index_i;
        // Walk up to the root and report where we are
        while (at_id != 0) {
            if (at_snarl) {
                at_id = dist_index.snarl_indexes[at_rank].parent_id;
                if (at_id != 0) {
                    at_snarl = !dist_index.snarl_indexes[at_rank].in_chain;
                    at_rank = at_snarl ? dist_index.get_primary_rank(at_id) : dist_index.get_chain_rank(at_id);
                }
            } else {
                at_id = dist_index.chain_indexes[at_rank].parent_id;
                if (at_id != 0) {
                    at_snarl = true;
                    at_rank = dist_index.get_primary_rank(at_id);
                }
            }
            if (at_id != 0) {
                if (at_snarl) {
                    cerr << "\tWhich is in parent snarl " << dist_index.snarl_indexes[at_rank].id_in_parent << "-" << dist_index.snarl_indexes[at_rank].end_id << endl;
                } else {
                    cerr << "\tWhich is in parent chain " << dist_index.chain_indexes[at_rank].id_in_parent << "-" << dist_index.chain_indexes[at_rank].end_id << endl;
                }
            }
        }
        cerr << "\tWhich is top level" << endl;
#endif

        //Keep track of all clusters on this snarl
        NodeClusters snarl_clusters(tree_state.all_seeds->size());

        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group, size_t& fragment_combined_group,
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
                    if (new_g != new_group) snarl_clusters.read_cluster_heads.erase(make_pair(read_num,new_group));
                    if (new_g != combined_group) snarl_clusters.read_cluster_heads.erase(make_pair(read_num,combined_group));
                    
                    snarl_clusters.read_cluster_heads.emplace(read_num,new_g);
                    end_dists = make_pair( min_not_minus_one(end_dists.first, old_dists.first),
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
        vector<pair<NetgraphNode, NodeClusters>>& child_nodes = tree_state.snarl_to_nodes[snarl_index_i];
        int64_t start_length = snarl_index.node_length(0);
        int64_t end_length = snarl_index.node_length(snarl_index.num_nodes*2 -1);


        //Maps each cluster of child nodes to its left and right distances
        //of the node its on
        hash_map<pair<size_t,size_t>, pair<int64_t, int64_t>> old_dists;
        old_dists.reserve(child_nodes.size());

        for (size_t i = 0; i < child_nodes.size() ; i++) {
            //Go through each child node of the netgraph

            NetgraphNode& child = child_nodes [i].first;

            // Get the node id of this netgraph node in its parent snarl
            // Ranks in parents are computed from node ID, so we have to get it.
            id_t child_node_id = child.id_in_parent(dist_index);

            //Rank of this node in the snarl
            //Note, if this node is a snarl/chain, then this snarl will be the secondary snarl
            size_t node_rank = child.rank_in_parent(dist_index, child_node_id);
            size_t rev_rank = node_rank % 2 == 0 ? node_rank + 1 : node_rank - 1;

            if (child.node_type == NODE) {
                //If this node is a node, we need to find the clusters
                int64_t node_len = snarl_index.node_length(node_rank);

                child_nodes[i].second = cluster_one_node(tree_state, child_node_id, node_len);

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
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
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

                pair<int64_t, int64_t> new_dists = snarl_index.dist_to_ends(node_rank,
                                        dists_c.first,dists_c.second);
#ifdef DEBUG_CLUSTER
cerr << "\tcluster: " << c_i << "dists to ends in snarl" << snarl_index.id_in_parent
     << " : " << new_dists.first << " " << new_dists.second << endl;
#endif

                snarl_clusters.fragment_best_left =min_not_minus_one( snarl_clusters.fragment_best_left,new_dists.first);
                snarl_clusters.fragment_best_right = min_not_minus_one(snarl_clusters.fragment_best_right, new_dists.second);
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
                size_t other_rank = other_node.rank_in_parent(dist_index, other_node_id);
                size_t other_rev = other_rank % 2 == 0 ? other_rank + 1 : other_rank - 1;

#ifdef DEBUG_CLUSTER
                cerr << "Other net graph node is " << typeToString(other_node.node_type)
                    << " headed by node " << other_node_id;


#endif


                //Find distance from each end of current node (i) to
                //each end of other node (j)
                int64_t dist_l_l = snarl_index.snarl_distance(rev_rank, other_rank);
                int64_t dist_l_r = snarl_index.snarl_distance(rev_rank, other_rev);
                int64_t dist_r_l = snarl_index.snarl_distance(node_rank, other_rank);
                int64_t dist_r_r = snarl_index.snarl_distance(node_rank, other_rev);

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
                         MinimumDistanceIndex::min_pos({dist_l_l, dist_l_r, dist_r_l, dist_r_r})-2 <= tree_state.read_distance_limit
                   && min_not_minus_one(curr_child_clusters.fragment_best_left, curr_child_clusters.fragment_best_right)-2
                                                <= tree_state.read_distance_limit) ||
                       (tree_state.fragment_distance_limit != 0 &&
                            MinimumDistanceIndex::min_pos({dist_l_l, dist_l_r,dist_r_l, dist_r_r})-2 <= tree_state.fragment_distance_limit
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


                        if (dist_l_l != -1 && dists_c.first != -1 && other_node_clusters.fragment_best_left != -1 ) {
                            //If cluster child_cluster_head can be combined with clusters in j
                            //from the left of both of them
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == -1 ? -1 :
                                  dist_l_l + dists_c.first + other_node_clusters.read_best_left[read_num]-1;
                            int64_t fragment_dist = dist_l_l + dists_c.first + other_node_clusters.fragment_best_left-1;
                            combine_clusters(c_group, group_l_l[read_num], fragment_group_l_l,
                                  fragment_dist, read_dist,  read_num);
                        }

                        if (dist_l_r != -1 && dists_c.first != -1 && other_node_clusters.fragment_best_right != -1 ) {
                            //If it can be combined from the left to the right of j
                            int64_t fragment_dist = dist_l_r + dists_c.first + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == -1 ? -1 :
                                 dist_l_r + dists_c.first + other_node_clusters.read_best_right[read_num]-1;
                            combine_clusters(c_group, group_l_r[read_num], fragment_group_l_r,
                                 fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != -1 && dists_c.second != -1 && other_node_clusters.fragment_best_left != -1 ) {
                            int64_t fragment_dist = dist_r_l + dists_c.second + other_node_clusters.fragment_best_left-1;
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == -1 ? -1 :
                                dist_r_l + dists_c.second + other_node_clusters.read_best_left[read_num]-1;
                            combine_clusters(c_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist,  read_num);
                        }
                        if (dist_r_r != -1 && dists_c.second != -1 && other_node_clusters.fragment_best_right != -1 ) {
                            int64_t fragment_dist = dist_r_r + dists_c.second + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == -1 ? -1 :
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


                        if (dist_l_l != -1 && curr_child_clusters.fragment_best_left != -1 && dists_k.first != -1 ){

                            int64_t fragment_dist = dist_l_l + curr_child_clusters.fragment_best_left + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == -1 ? -1 :
                                dist_l_l + curr_child_clusters.read_best_left[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_l_l[read_num], fragment_group_l_l,
                                fragment_dist,read_dist, read_num);
                        }
                        if (dist_l_r != -1 && curr_child_clusters.fragment_best_left != -1 && dists_k.second != -1  ) {

                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == -1 ? -1 :
                               dist_l_r + curr_child_clusters.read_best_left[read_num] + dists_k.second-1;
                            int64_t fragment_dist = dist_l_r + curr_child_clusters.fragment_best_left + dists_k.second-1;
                            combine_clusters(k_group, group_l_r[read_num], fragment_group_l_r,
                               fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != -1 && curr_child_clusters.fragment_best_right != -1 && dists_k.first != -1  ) {

                            int64_t fragment_dist = dist_r_l + curr_child_clusters.fragment_best_right + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == -1 ? -1 :
                                dist_r_l + curr_child_clusters.read_best_right[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_r != -1 && curr_child_clusters.fragment_best_right != -1 && dists_k.second != -1 ) {

                            int64_t fragment_dist = dist_r_r + curr_child_clusters.fragment_best_right + dists_k.second-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == -1 ? -1 :
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



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_chain(TreeState& tree_state, size_t chain_i, size_t depth) const {

        //Get the index of the chain we're on
        MinimumDistanceIndex::ChainIndex& chain_index = dist_index.chain_indexes[chain_i];

#ifdef DEBUG_CLUSTER
        cerr << "Cluster chain " << chain_index.id_in_parent << endl;
#endif

        //Maps each snarl to its clusters, in the order of the snarls in the chain
        //TODO: I think this doesn't need to be a map anymore since we're sorting later anyway, depending on the sort algorithm
        hash_map<size_t, pair<size_t, NodeClusters>>& snarls_in_chain =tree_state.chain_to_snarls[chain_i];

        size_t connected_component_num = dist_index.get_connected_component(chain_index.id_in_parent);
        int64_t chain_length = depth == 0 ? dist_index.top_level_chain_length(chain_index.id_in_parent)
                                          : chain_index.chain_length();

        /*Start by making a list of all clusters (snarl and seed) ordered by their occurrence in the chain
         * Snarl clusters are only ordered relative to seed clusters, and we never compare snarl clusters
         * that came from the same snarl cluster
         * When we add the snarl clusters to cluster_head_indices, check if clusters in the same snarl can combine
         * and update each one's distance to the ends of the chain
         */

        //Holds all clusters as tuple of <is top level seed ? -1 : offset of first node in the snarl, read num, seed num> 
        //vectors of tuples are sorted by the offset
        vector<tuple<int64_t, size_t, size_t>> cluster_head_indices ;
        cluster_head_indices.reserve(snarls_in_chain.size());


        //Insert a cluster (offset in chain, read num, seed num) into a sorted list of clusters
        auto insert_in_order = [&](vector<tuple<int64_t, size_t, size_t>>& input_vector, tuple<int64_t, size_t, size_t> item) {
            //Get the offset of the thing we're inserting
            const Seed& curr_seed = tree_state.all_seeds->at(std::get<1>(item))->at(std::get<2>(item));
            pos_t pos =  curr_seed.pos;
            //Offset in the chain (if its a top-level seed) or the offset of the end of the first boundary node in the chain
            //of a snarl if its a snarl cluster
            //snarls always go after seeds on the same node since seeds will be added first
            int64_t item_offset = std::get<0>(item) == -1 ? curr_seed.offset : std::get<0>(item);

            //Get an iterator to where the thing we're inserting should go in the vector
            std::vector<tuple<int64_t, size_t, size_t>>::iterator insert_itr = std::upper_bound(input_vector.begin(), input_vector.end(),
               item_offset, [&](int64_t val, tuple<int64_t, size_t, size_t> vector_item) {
                   int64_t offset = std::get<0>(vector_item) == -1 ? 
                        tree_state.all_seeds->at(std::get<1>(vector_item))->at(std::get<2>(vector_item)).offset :
                        std::get<0>(vector_item);
                   return val < offset; 
               });
            //Add the item to the vector in sorted order
            input_vector.insert(insert_itr, item);
        };

        //Add top-level seed clusters and top-level simple snarl clusters to sorted list
        if (depth == 0 && tree_state.component_to_index.count(connected_component_num) != 0) {

            //Add the top-level seeds
            for (const pair<size_t, size_t>& indices : tree_state.top_level_seed_clusters[tree_state.component_to_index[connected_component_num]]) {
                insert_in_order(cluster_head_indices, make_tuple(-1, indices.first, indices.second));
            }

            //Cluster top-level simple snarls and add them to the list of cluster heads after the seeds
            for (auto& snarl_to_node : 
                    tree_state.simple_snarl_to_nodes_by_component[tree_state.component_to_index[connected_component_num]] ){

                size_t start_rank = snarl_to_node.first;
                bool rev_in_chain; id_t node_id; int64_t node_len; //Don't really need these
                int64_t start_length ; int64_t end_length;
                std::tie(node_id, rev_in_chain, start_length, end_length, node_len) = snarl_to_node.second[0];
                int64_t snarl_length = start_rank == 0 ? chain_index.prefix_sum[start_rank + 1] - start_length : 
                        chain_index.prefix_sum[start_rank + 1] - chain_index.prefix_sum[start_rank] - start_length; 

                //Loop distances from and to the ends of the nodes in the snarl
                int64_t loop_left = chain_index.loop_rev[start_rank] == 0 ? -1 : chain_index.loop_rev[start_rank] - 1 + start_length;
                int64_t loop_right = chain_index.loop_fd[start_rank+1] == 0 ? -1 : chain_index.loop_fd[start_rank+1] - 1 + end_length;
                if (rev_in_chain) {
                    int64_t tmp = loop_left;
                    loop_left = loop_right;
                    loop_right = tmp;
                }
                
                // Grab the snarl's overall rank so we can identify it.
                size_t snarl_rank = dist_index.get_primary_assignment(node_id);

                //Cluster this snarl
                //Updates distances to the ends of the node (not including boundary nodes of the snarl), sides
                //are relative to the node, not orientation in the chain
                hash_set<pair<size_t, size_t>> simple_snarl_clusters = 
                        cluster_simple_snarl(tree_state, snarl_rank, snarl_to_node.second, loop_left, loop_right, snarl_length); 

                //The offset of this snarl in the chain for placement in cluster_head_indices
                int64_t offset = start_rank == 0 ? start_length : chain_index.prefix_sum[start_rank] + start_length - 1;

                int64_t add_dist_left_left = start_rank == 0 ? start_length : chain_index.prefix_sum[start_rank] - 1 + start_length;
    
                int64_t add_dist_right_right = chain_length+1 -chain_index.prefix_sum[start_rank+1];
#ifdef DEBUG_CLUSTER
                cerr << "\tupdating simple snarl distances to the ends of the chain with additional distances " 
                     << add_dist_left_left << " " << add_dist_right_right << endl;
                MinimumDistanceIndex::SnarlIndex& snarl_index = dist_index.snarl_indexes[snarl_rank];
                assert(dist_index.get_chain_rank(snarl_index.id_in_parent) ==  start_rank);
                assert(snarl_index.node_length(dist_index.get_primary_rank(node_id)) == node_len);
                if (snarl_index.rev_in_parent) {
                    assert(snarl_index.node_length(0) == end_length);
                    assert(snarl_index.node_length(snarl_index.num_nodes*2-1) == start_length);
                } else {
                    assert(snarl_index.node_length(0) == start_length);
                    assert(snarl_index.node_length(snarl_index.num_nodes*2-1) == end_length);
                }

#endif
    
                for (const pair<size_t, size_t>& cluster_head : simple_snarl_clusters) {
                    //Add this new snarl cluster to list. It will be treated as any other snarl cluster

                    insert_in_order(cluster_head_indices, make_tuple( offset, cluster_head.first, cluster_head.second));

                    int64_t old_left = rev_in_chain ? tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second : 
                                                        tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first;
                    int64_t old_right = rev_in_chain ? tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first : 
                                                        tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second;

    
                    //Update distances to reach the ends of the chain
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first  = 
                        old_left == -1 ? -1 : old_left +  add_dist_left_left; 
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = 
                        old_right == -1 ? -1 : old_right + add_dist_right_right; 
#ifdef DEBUG_CLUSTER
                        cerr << "\tkeeping simple snarl cluster: " 
                            << tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).pos << " with distances: " 
                            << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first << ", " 
                            << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second << endl;
#endif

                }
            }
        }

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
                          min_not_minus_one(old_dists.first, dists.first),
                          min_not_minus_one(old_dists.second, dists.second));
                    tree_state.read_cluster_dists[read_num][combined_group] = dists;
#ifdef DEBUG_CLUSTER
                    cerr << " New dists for read num " << read_num << ": "
                         << tree_state.read_cluster_dists[read_num][combined_group].first << " "
                         << tree_state.read_cluster_dists[read_num][combined_group].second << endl;
#endif
                }

                if (tree_state.fragment_distance_limit != 0 && fragment_dist != -1) {
                    if (fragment_combined_group != -1) {
                    //If we're also keeping track of fragment clusters
                        tree_state.fragment_union_find.union_groups(fragment_combined_group,
                                                            new_group + tree_state.read_index_offsets[read_num]);
                    }
                    fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                         new_group + tree_state.read_index_offsets[read_num]);
                }
            } else if (tree_state.fragment_distance_limit != 0 && fragment_dist != -1 &&
                        fragment_dist <= tree_state.fragment_distance_limit) {
                //If these aren't in the same read cluster but are in
                //the same fragment cluster
                if (fragment_combined_group != -1) {
                    tree_state.fragment_union_find.union_groups(
                                    fragment_combined_group, new_group + tree_state.read_index_offsets[read_num]);
                }
                fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                     new_group + tree_state.read_index_offsets[read_num]);
            }
            return;
        };


        //Go through the chain by snarl. 
        //Potentially clusters clusters in the same snarl if there is a path between them in the chain
        //update the distances to the ends of the chain
        for (auto& kv: snarls_in_chain) {
            /*
             * Snarls are in the order that they are traversed in the chain
             * snarls_in_chain is a map from rank of a snarl to the snarl index
             *  in dist_index.snarl_indexes and NodeClusters for that snarl
             */


            //rank of the boundary node of the snarl that occurs first in
            //the chain
            size_t start_rank = kv.first;

            //The clusters of the current snarl
            NodeClusters& snarl_clusters = kv.second.second;

            //Index of the current snarl
            MinimumDistanceIndex::SnarlIndex& snarl_index = dist_index.snarl_indexes[kv.second.first];

            int64_t start_length = snarl_index.rev_in_parent
                         ? snarl_index.node_length(snarl_index.num_nodes * 2 - 1)
                         : snarl_index.node_length(0);
            int64_t end_length = snarl_index.rev_in_parent ? snarl_index.node_length(0)
                        : snarl_index.node_length(snarl_index.num_nodes * 2 - 1);
            int64_t snarl_length = snarl_index.snarl_length();


            //Combine snarl clusters that can be reached by looping
            int64_t loop_dist_end = chain_index.loop_fd[start_rank + 1] - 1 ;
            int64_t loop_dist_start = chain_index.loop_rev[start_rank] - 1;


            //Distance from the ends of the current snarl to the ends of the chain
            int64_t add_dist_left_left = start_rank == 0 ? 0 : chain_index.prefix_sum[start_rank] - 1;

            int64_t add_dist_right_right = start_rank + 1 == chain_index.prefix_sum.size() - 2 ? 0 : 
                                    chain_index.prefix_sum[chain_index.prefix_sum.size()-1] - chain_index.prefix_sum[start_rank+1] - end_length;

            int64_t add_dist_left_right = start_rank == 0 || loop_dist_start == -1 ? -1 : 
                                          loop_dist_start + add_dist_right_right + snarl_length - start_length;
            int64_t add_dist_right_left = start_rank+1 == chain_index.prefix_sum.size() - 2 || loop_dist_end == -1 ? -1 :
                                           loop_dist_end+ add_dist_left_left + snarl_length - end_length;

            hash_set<pair<size_t,size_t>> to_add;//new cluster group ids from snarl clusters
            vector<pair<size_t,size_t>> to_erase; //old cluster group ids

            //Combined snarl clusters by taking chain loop left/right
            vector<size_t> snarl_cluster_left (tree_state.all_seeds->size(),-1);
            vector<size_t> snarl_cluster_right (tree_state.all_seeds->size(), -1);
            size_t fragment_snarl_cluster_left = -1;
            size_t fragment_snarl_cluster_right = -1;

            for (pair<size_t, size_t> cluster_head : snarl_clusters.read_cluster_heads) {
                // For each of the clusters for the current snarl,
                // first check if it can be combined with another cluster
                // in the same snarl by taking loops in the chain,
                // then, find if it belongs to the new combined cluster
                // that includes chain clusters
                size_t read_num = cluster_head.first;

                pair<int64_t, int64_t> snarl_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if (loop_dist_start != -1) {
                    //If there is a loop going out and back into the start of
                    //the snarl, this cluster may be combined with other snarl clusters

                    //The distance to the right side of the snarl that is found by taking the leftmost seed and
                    // looping through the chain to the left
                    int64_t new_right = snarl_dists.first == -1 || loop_dist_start == -1
                                        ? -1
                                        : snarl_dists.first + loop_dist_start + snarl_length - start_length;
                    if (snarl_dists.second == -1 || (new_right != -1 && new_right < snarl_dists.second)){
                        snarl_dists.second = min_not_minus_one(new_right, snarl_dists.second);
                        snarl_clusters.fragment_best_right = min_not_minus_one(snarl_clusters.fragment_best_right, new_right);
                        snarl_clusters.read_best_right[read_num] = min_not_minus_one(snarl_clusters.read_best_right[read_num], new_right);
#ifdef DEBUG_CLUSTER
cerr << "  Updating looping distance to right of snarl cluster " << read_num <<":" << cluster_head.second << ": "
     << new_right << " -> " << snarl_dists.second <<  endl;
#endif
                    }


                    if (snarl_clusters.fragment_best_left!= -1 && snarl_dists.first != -1 ) {
                        //If this cluster can be combined with another cluster
                        //from the left

#ifdef DEBUG_CLUSTER
cerr << "  Combining this cluster from the left " << endl;
#endif
                        int64_t read_dist =  snarl_clusters.read_best_left[read_num] == -1 ? -1 :  
                                    snarl_clusters.read_best_left[read_num] + snarl_dists.first + loop_dist_start - start_length - 1;
                        int64_t fragment_dist = snarl_clusters.fragment_best_left == -1 ? -1 :
                                    snarl_clusters.fragment_best_left + snarl_dists.first + loop_dist_start - start_length - 1;

                        combine_snarl_clusters(cluster_head.second, snarl_cluster_left[read_num], 
                                fragment_snarl_cluster_left,  to_erase, fragment_dist, read_dist, snarl_dists, read_num);
                    }

                }

                if (loop_dist_end != -1) {
                    //If there is a loop to the right
                    int64_t new_left = snarl_dists.second == -1 || loop_dist_end == -1 ? -1
                          : snarl_dists.second + loop_dist_end + snarl_length - end_length;
                    if (snarl_dists.first == -1 || (new_left != -1 && new_left < snarl_dists.first)){
                        //If this is an improvement, update distances
                        snarl_dists.first = new_left;
                        snarl_clusters.read_best_left[read_num] = min_not_minus_one(new_left, snarl_clusters.read_best_left[read_num]);
                        snarl_clusters.fragment_best_left = min_not_minus_one(new_left, snarl_clusters.fragment_best_left);

#ifdef DEBUG_CLUSTER
cerr << "Updating looping distance to left of snarl cluster " << read_num << ":" << cluster_head.second << ": "
     << new_left << endl;
#endif
                    }

                    if (snarl_clusters.fragment_best_right != -1 && snarl_dists.second != -1) {
                        //If this cluster can be combined with another cluster
                        //from the right

#ifdef DEBUG_CLUSTER
cerr << "  Maybe combining this cluster from the right" << endl;
#endif
                        int64_t read_dist = snarl_clusters.read_best_right[read_num] == -1 ? -1 :
                            snarl_clusters.read_best_right[read_num] + snarl_dists.second  + loop_dist_end - end_length - 1;
                        int64_t fragment_dist = snarl_clusters.fragment_best_right == -1 ? -1 : 
                            snarl_clusters.fragment_best_right + snarl_dists.second + loop_dist_end - end_length - 1;

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
            int64_t offset_in_chain = start_rank == 0 ? start_length : chain_index.prefix_sum[start_rank] + start_length - 1;
            for (auto& cluster_head : to_add) {
                //Add the clusters on this snarl to our overall list of clusters and update the distances
                //for each of the clusters
                insert_in_order(cluster_head_indices, make_tuple(offset_in_chain, cluster_head.first, cluster_head.second));

                //Get the distance to the start side of the chain
                int64_t dist_left_left = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first == -1 ? -1 
                        : tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first +  add_dist_left_left;
                int64_t dist_right_left = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second == -1 
                        || add_dist_right_left == -1 ? -1 
                        : tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second +  add_dist_right_left;


                int64_t dist_right_right = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second == -1 ? -1 
                        : tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second + add_dist_right_right;
                int64_t dist_left_right = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first == -1 ||
                        add_dist_left_right == -1 ? -1 
                        : tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first + add_dist_left_right;

                tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first  = min_not_minus_one(dist_left_left, dist_right_left); 
                tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = min_not_minus_one(dist_right_right, dist_left_right); 

            }

        }



        //Keep track the snarl clusters from the current snarl and all previous snarl clusters
        //since the last time we saw a top-level seed
        //(We don't necessarily want to combine two clusters if they aren't top level seeds, since they would have
        //already been clustered if they were reachable, whereas top-level seeds are always reachable)
        //cluster head, distance left, distance right
        vector<tuple<size_t, int64_t, int64_t>> prev_snarl_cluster_fragment;

        vector<tuple<size_t, int64_t, int64_t>> snarl_cluster_fragment;

        //The rightmost seed that is a top-level seed or in a cluster with a top-level seed
        //<cluster head, distance to the left of the chain
        pair<size_t, int64_t> seed_cluster_fragment (0, -1);


        //And the same thing for each read
        vector<vector<tuple<size_t, int64_t, int64_t>>> snarl_cluster_by_read(tree_state.all_seeds->size());

        vector<vector<tuple<size_t, int64_t, int64_t>>> prev_snarl_cluster_by_read(tree_state.all_seeds->size());

        vector<pair<size_t, int64_t>> seed_cluster_by_read (tree_state.all_seeds->size(), make_pair(0, -1));

        //The offset of the last snarl we got clusters from
        vector<int64_t> last_snarl_rank (tree_state.all_seeds->size(), -1);
        int64_t last_snarl_rank_fragment = -1;


        //The clusters of the chain that are built from the snarl and minimizer clusters
        //This will get updated as we traverse through the child clusters
        NodeClusters chain_clusters(tree_state.all_seeds->size());

        //Go through the clusters in order and cluster them. If the cluster was a real cluster, then compare it to the
        //last seed cluster we found. If it was a seed, compare it to all previous clusters since we last saw a seed
        for (tuple<int64_t, size_t, size_t>& seed_index : cluster_head_indices) {

            vector<pair<size_t, size_t>> to_erase;

            //The new cluster head for whatever this ends up clustering with
            size_t new_cluster_head = std::get<2>(seed_index);

            size_t read_num = std::get<1>(seed_index);
            if (std::get<0>(seed_index) == -1) {

                //If this is a top-level seed, then try to compare it to the most recent snarl clusters and update the
                //best seed-only clusters. Also compare it to the last seed cluster

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

                chain_clusters.read_best_left[read_num] = min_not_minus_one(chain_clusters.read_best_left[read_num], offset);
                chain_clusters.read_best_right[read_num] = min_not_minus_one(chain_clusters.read_best_right[read_num], right_offset);
                chain_clusters.fragment_best_left = min_not_minus_one( chain_clusters.fragment_best_left, offset);
                chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right, right_offset);

                last_snarl_rank_fragment = -1;
                last_snarl_rank[read_num] = -1;


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
                            if (std::get<2>(prev_snarl_cluster_fragment[i]) != -1 && right_offset != -1 && distance_in_chain <= tree_state.fragment_distance_limit) {
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
#ifdef DEBUG_CLUSTER
                            else if (std::get<2>(prev_snarl_cluster_fragment[i]) != -1 && right_offset != -1) {
                                cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << std::get<0>(prev_snarl_cluster_fragment[i]) << " is too large" << endl;
                            }
#endif
                        }
                    }

                    //Combine this seed cluster with the last seed cluster we found
#ifdef DEBUG_CLUSTER
                    cerr << "\t Compare to last seed cluster fragment: " << endl;
#endif
                    int64_t distance_in_chain = offset - seed_cluster_fragment.second;
#ifdef DEBUG_CLUSTER
                    //assert(distance_in_chain >= 0);
#endif
                    if (seed_cluster_fragment.second != -1 && offset != -1 && distance_in_chain <= tree_state.fragment_distance_limit) {
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
#ifdef DEBUG_CLUSTER
                    else if (seed_cluster_fragment.second != -1 && offset != -1) {
                        cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << seed_cluster_fragment.first << " is too large" << endl;
                    }
#endif
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
                        cerr << "\t\t cluster" << read_num << " " << std::get<0>(prev_snarl_cluster_by_read[read_num][i]) << " " << 
                                tree_state.all_seeds->at(read_num)->at(std::get<0>(prev_snarl_cluster_by_read[read_num][i])).pos << 
                                "(distances: " <<  std::get<1>(prev_snarl_cluster_by_read[read_num][i]) << ", " << std::get<2>(prev_snarl_cluster_by_read[read_num][i]) << ")";
#endif
                        if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != -1 && right_offset != -1 && 
                                distance_in_chain <= tree_state.read_distance_limit) {
                            
                            tree_state.read_union_find[read_num].union_groups(std::get<2>(seed_index), 
                                                                              std::get<0>(prev_snarl_cluster_by_read[read_num][i])); 
                            size_t new_group = tree_state.read_union_find[read_num].find_group(std::get<2>(seed_index));
                            if (new_group == new_cluster_head) {
                                to_erase.emplace_back(read_num, std::get<0>(prev_snarl_cluster_by_read[read_num][i]));
                            } else {
                                to_erase.emplace_back(read_num, new_cluster_head);
                            }
                            best_left = min_not_minus_one(best_left, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].first);
                            best_right = min_not_minus_one(best_right, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].second);
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
#ifdef DEBUG_CLUSTER
                        else if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != -1 && right_offset != -1) {
                            cerr << "\t\tDistance " << distance_in_chain << " to it is too large" << endl;
                        }
#endif
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
                if (seed_cluster_by_read[read_num].second != -1 && offset != -1 && 
                        distance_in_chain <= tree_state.read_distance_limit) {
                    tree_state.read_union_find[read_num].union_groups(std::get<2>(seed_index), 
                                                                          seed_cluster_by_read[read_num].first); 
                    size_t new_group = tree_state.read_union_find[read_num].find_group(std::get<2>(seed_index));
                    if (new_group == new_cluster_head) {
                        to_erase.emplace_back(read_num, seed_cluster_by_read[read_num].first);
                    } else {
                        to_erase.emplace_back(read_num, new_cluster_head);
                    }
                    best_left = min_not_minus_one(best_left, seed_cluster_by_read[read_num].second);
                    best_right = min_not_minus_one(best_right, chain_length-seed_cluster_by_read[read_num].second+1);
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
#ifdef DEBUG_CLUSTER
                else if (seed_cluster_by_read[read_num].second != -1 && offset != -1) {
                    cerr << "\t\tDistance " << distance_in_chain << " to it is too large" << endl;
                }
#endif

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

                chain_clusters.read_best_left[read_num] = min_not_minus_one(chain_clusters.read_best_left[read_num], offset);
                chain_clusters.read_best_right[read_num] = min_not_minus_one(chain_clusters.read_best_right[read_num], right_offset);
                chain_clusters.fragment_best_left = min_not_minus_one( chain_clusters.fragment_best_left, offset);
                chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right, right_offset);
#ifdef DEBUG_CLUSTER
                cerr << "At a snarl cluster on read " << read_num << " with head " << tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                        " and distances " << offset << " " << right_offset << endl;
#endif

                if (std::get<0>(seed_index) != last_snarl_rank[read_num]) {

                    //If we're moving on to a new snarl, update the most recent snarl clusters
                    prev_snarl_cluster_by_read[read_num].insert(prev_snarl_cluster_by_read[read_num].end(),
                        snarl_cluster_by_read[read_num].begin(), snarl_cluster_by_read[read_num].end());

                    last_snarl_rank[read_num] = std::get<0>(seed_index);
                }


                //try clustering by fragment
                if (tree_state.fragment_distance_limit != 0) {
                    int64_t new_best_right = -1;
                    if (std::get<0>(seed_index) != last_snarl_rank_fragment) {
#ifdef DEBUG_CLUSTER
                        cerr << "\tOn a new snarl with offset in chain " << std::get<0>(seed_index) << endl;
#endif
                        //And the same for the fragment clusters
                        prev_snarl_cluster_fragment.insert(prev_snarl_cluster_fragment.end(),
                            snarl_cluster_fragment.begin(), snarl_cluster_fragment.end());
                        snarl_cluster_fragment.clear();
    
    
                        last_snarl_rank_fragment = std::get<0>(seed_index);
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
                            if (std::get<2>(prev_snarl_cluster_fragment[i])  != -1 && offset != -1 && 
                                    distance_in_chain <= tree_state.fragment_distance_limit) {
                                tree_state.fragment_union_find.union_groups(std::get<0>(prev_snarl_cluster_fragment[i]), 
                                                        std::get<2>(seed_index)+tree_state.read_index_offsets[read_num]);
                                new_best_right = min_not_minus_one(new_best_right, std::get<2>(prev_snarl_cluster_fragment[i]));
#ifdef DEBUG_CLUSTER
                                cerr << "\t\tCombining fragment: snarl cluster on component " << connected_component_num << ", " << 
                                        tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                        "(distances: " << offset << ", " << right_offset << ")" <<  
                                        " with snarl cluster" << " " << read_num << " " << std::get<0>(prev_snarl_cluster_fragment[i]) << " " << 
                                        "(distances: " << std::get<1>(prev_snarl_cluster_fragment[i]) << ", " <<  std::get<2>(prev_snarl_cluster_fragment[i]) << 
                                        ") with distance between them " << distance_in_chain << endl;
#endif
                            } 
#ifdef DEBUG_CLUSTER
                            else if(std::get<2>(prev_snarl_cluster_fragment[i])  != -1 && offset != -1) {
                                cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << std::get<0>(prev_snarl_cluster_fragment[i]) << " is too large" << endl;
                            }
#endif
                        }
                    } 
                    //Always try to compare this to the most recent seed cluster
#ifdef DEBUG_CLUSTER
                    cerr << "\tCompare to last seed cluster fragment:" << endl;
#endif
                    int64_t distance_in_chain = offset - seed_cluster_fragment.second;
                    if (seed_cluster_fragment.second != -1 && offset != -1 && distance_in_chain <= tree_state.fragment_distance_limit) {
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
#ifdef DEBUG_CLUSTER
                    else if (seed_cluster_fragment.second != -1 && offset != -1) {
                        cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << seed_cluster_fragment.first << " is too large" << endl;
                    }
#endif
                    snarl_cluster_fragment.emplace_back(std::get<2>(seed_index)+tree_state.read_index_offsets[read_num],
                        offset, min_not_minus_one(new_best_right, right_offset));
                }


                //try clustering by read
                if (!prev_snarl_cluster_by_read[read_num].empty()) {
#ifdef DEBUG_CLUSTER
                    cerr << "\tComparing to last snarl clusters read "  << endl;
#endif
                    //If we saw a snarl cluster in a previous snarl and haven't seen a seed cluster since then
                    for (size_t i = 0 ; i < prev_snarl_cluster_by_read[read_num].size() ; i++) {
                        int64_t distance_in_chain = offset + std::get<2>(prev_snarl_cluster_by_read[read_num][i]) - chain_length - 1;
                        if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != -1 && offset != -1 &&
                                distance_in_chain <= tree_state.read_distance_limit) {
                            tree_state.read_union_find[read_num].union_groups(std::get<0>(prev_snarl_cluster_by_read[read_num][i]), 
                                                                          new_cluster_head);
#ifdef DEBUG_CLUSTER
                            cerr << "\t\tSnarl cluster " << read_num << " " << std::get<0>(prev_snarl_cluster_by_read[read_num][i]) << " " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<0>(prev_snarl_cluster_by_read[read_num][i])).pos << 
                                    " (distances: " << std::get<1>(prev_snarl_cluster_by_read[read_num][i]) << ", " <<  
                                     std::get<2>(prev_snarl_cluster_by_read[read_num][i]) <<  ")" << endl; 
#endif
                            size_t new_group = tree_state.read_union_find[read_num].find_group(new_cluster_head);
                            if (new_group == new_cluster_head) {
                                to_erase.emplace_back(read_num, std::get<0>(prev_snarl_cluster_by_read[read_num][i]));
                            } else {
                                to_erase.emplace_back(read_num, new_cluster_head);
                            }
                            best_left = min_not_minus_one(best_left, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].first);
                            best_right = min_not_minus_one(best_right, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].second);
                            tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right);
                            std::get<0>(prev_snarl_cluster_by_read[read_num][i]) = new_group;
                            new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                            cerr << "\t\t\t...Combining read: snarl cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(new_cluster_head).pos << 
                                    " (distances: " << offset << ", " << right_offset << ")" <<  
                                    " with snarl cluster with distance between them " << distance_in_chain<< endl;
#endif
                        } 
#ifdef DEBUG_CLUSTER
                        else if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != -1 && offset != -1) {
                            cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << std::get<0>(prev_snarl_cluster_by_read[read_num][i]) << " is too large" << endl;
                        }
#endif
                    }

                }
                //Always try to cluster with the last seed cluster
                int64_t distance_in_chain = offset - seed_cluster_by_read[read_num].second; 
#ifdef DEBUG_CLUSTER
                cerr << "\tComparing to last seed cluster read "  << read_num << " " << seed_cluster_by_read[read_num].first << " " <<  
                         tree_state.all_seeds->at(read_num)->at(seed_cluster_by_read[read_num].first).pos << 
                         " (distances: " << seed_cluster_by_read[read_num].second << ", " << chain_length - seed_cluster_by_read[read_num].second + 1 << ")" << endl;
#endif
                if (seed_cluster_by_read[read_num].second != -1 && offset != -1 && 
                        distance_in_chain <= tree_state.read_distance_limit) {

                    tree_state.read_union_find[read_num].union_groups(seed_cluster_by_read[read_num].first, 
                                                                      new_cluster_head);
                    size_t new_group = tree_state.read_union_find[read_num].find_group(new_cluster_head);
                    if (new_group == new_cluster_head) {
                        to_erase.emplace_back(read_num, seed_cluster_by_read[read_num].first);
                    } else {
                        to_erase.emplace_back(read_num, new_cluster_head);
                    }
                    best_left = min_not_minus_one(best_left, seed_cluster_by_read[read_num].second);
                    best_right = min_not_minus_one(best_right, chain_length-seed_cluster_by_read[read_num].second+1);
                    tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right);
                    seed_cluster_by_read[read_num].first = new_group;
                    new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                        cerr << "\t\tCombining read: snarl cluster on component " << connected_component_num << ", " << 
                                tree_state.all_seeds->at(read_num)->at(new_cluster_head).pos << 
                                "(distances: " << offset << ", " << right_offset << ")" <<  
                                " with seed cluster with distance between them " << distance_in_chain << endl;
#endif

                    right_offset = min_not_minus_one(right_offset, chain_length-seed_cluster_by_read[read_num].second+1);
                } 
#ifdef DEBUG_CLUSTER
                else if (seed_cluster_by_read[read_num].second != -1 && offset != -1) {
                    cerr << "\t\tDistance " << distance_in_chain << " to " << read_num << " " << seed_cluster_by_read[read_num].first << " is too large" << endl;
                }
#endif

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
            vector<size_t> combined_cluster (tree_state.all_seeds->size(), -1);
            size_t fragment_combined_cluster = -1;

            for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
                //For each chain cluster
                size_t read_num = cluster_head.first;
                pair<int64_t, int64_t>& chain_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if ((chain_dists.second != -1 && chain_clusters.read_best_left[read_num] != -1 &&
                     chain_dists.second + chain_clusters.read_best_left[read_num] - first_length - 1 <= tree_state.read_distance_limit) ||
                   (chain_dists.first != -1 && chain_clusters.read_best_right[read_num] != -1 &&
                      chain_dists.first + chain_clusters.read_best_right[read_num] - first_length - 1 <= tree_state.read_distance_limit)){
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
                        combined_cluster[read_num] = new_group;
                    }

                    if (tree_state.fragment_distance_limit != 0) {
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second + tree_state.read_index_offsets[read_num]);
                    }
                } else if (tree_state.fragment_distance_limit != 0 &&
                   ((chain_dists.second != -1 && chain_clusters.fragment_best_left != -1 &&
                     chain_dists.second + chain_clusters.fragment_best_left - first_length - 1 <= tree_state.fragment_distance_limit) ||
                   (chain_dists.first != -1 && chain_clusters.fragment_best_right != -1 &&
                      chain_dists.first + chain_clusters.fragment_best_right - first_length - 1 <= tree_state.fragment_distance_limit))){
                    //If we can cluster by fragment
                    if (fragment_combined_cluster != -1) {
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
                        assert(dists.first == -1 || dists.first >= chain_clusters.read_best_left[read_num]);
                        assert(dists.second == -1 || dists.second >= chain_clusters.read_best_right[read_num]);
                        assert(dists.first == -1 || dists.first >= chain_clusters.fragment_best_left);
                        assert(dists.second == -1 || dists.second >= chain_clusters.fragment_best_right);
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
                assert(!any_clusters || got_read_left || chain_clusters.read_best_left[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_left[read_num] == -1);
                assert(!any_clusters || got_read_right || chain_clusters.read_best_right[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_right[read_num] == -1);
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


    hash_set<pair<size_t, size_t>> SnarlSeedClusterer::cluster_simple_snarl(TreeState& tree_state, size_t snarl_rank, vector<tuple<id_t, bool, int64_t, int64_t, int64_t>> nodes, 
                                int64_t loop_left, int64_t loop_right, int64_t snarl_length) const {
        //Cluster a top-level simple snarl and save the distances to the ends of the node in tree_state.read_cluster_dists 
        //Returns a vector of cluster heads (<read_num, seed num>)
        
        auto& snarl_index = dist_index.snarl_indexes[snarl_rank];
        
#ifdef DEBUG_CLUSTER
        cerr << "Cluster simple snarl " << snarl_index.id_in_parent << "-" << snarl_index.end_id << " which has rank " << snarl_rank << endl;
        cerr << " loop distances: " << loop_left << " " << loop_right  << " snarl length " << snarl_length << endl;
        
        bool at_snarl = true;
        id_t at_id = snarl_index.id_in_parent;
        size_t at_rank = snarl_rank;
        // Walk up to the root and report where we are
        while (at_id != 0) {
            if (at_snarl) {
                at_id = dist_index.snarl_indexes[at_rank].parent_id;
                if (at_id != 0) {
                    at_snarl = !dist_index.snarl_indexes[at_rank].in_chain;
                    at_rank = at_snarl ? dist_index.get_primary_rank(at_id) : dist_index.get_chain_rank(at_id);
                }
            } else {
                at_id = dist_index.chain_indexes[at_rank].parent_id;
                if (at_id != 0) {
                    at_snarl = true;
                    at_rank = dist_index.get_primary_rank(at_id);
                }
            }
            if (at_id != 0) {
                if (at_snarl) {
                    cerr << "\tWhich is in parent snarl " << dist_index.snarl_indexes[at_rank].id_in_parent << "-" << dist_index.snarl_indexes[at_rank].end_id << endl;
                } else {
                    cerr << "\tWhich is in parent chain " << dist_index.chain_indexes[at_rank].id_in_parent << "-" << dist_index.chain_indexes[at_rank].end_id << endl;
                }
            }
        }
        cerr << "\tWhich is top level" << endl;
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
            fragment_best_left = min_not_minus_one(fragment_best_left, clusters.fragment_best_left);
            fragment_best_right = min_not_minus_one(fragment_best_right, clusters.fragment_best_right);
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num ++) {
                best_left[read_num] = min_not_minus_one(best_left[read_num], clusters.read_best_left[read_num]);
                best_right[read_num] = min_not_minus_one(best_right[read_num], clusters.read_best_right[read_num]); 
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
                                min_not_minus_one(dist_right, dist_left + loop_left + snarl_length);
                    } else {
                        int64_t best_left = min_not_minus_one(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].first);
                        int64_t best_right =min_not_minus_one(dist_right, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].second);

                        tree_state.read_union_find[cluster_head.first].union_groups(cluster_head.second, combined_read_left[cluster_head.first]);
                        if (tree_state.read_union_find[cluster_head.first].find_group(cluster_head.second) == cluster_head.second) {
                            result.erase(make_pair(cluster_head.first, combined_read_left[cluster_head.first]));
                            combined_read_left[cluster_head.first] = cluster_head.second;
                        } else {
                            result.erase(cluster_head);
                        }


                        tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].first = 
                            min_not_minus_one(best_left, loop_right == -1 ? -1 : best_right + loop_right + snarl_length);
                        tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].second = 
                            min_not_minus_one(best_right, best_left + loop_left + snarl_length);
                    }
                } else {
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = 
                            min_not_minus_one(dist_right, dist_left + loop_left + snarl_length);
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
                                min_not_minus_one(dist_left, dist_right + loop_right + snarl_length);
                    } else {
                        int64_t best_left = min_not_minus_one(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].first);
                        int64_t best_right = min_not_minus_one(dist_right, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].second);
                        tree_state.read_union_find[cluster_head.first].union_groups(cluster_head.second, combined_read_right[cluster_head.first]);
                        if (tree_state.read_union_find[cluster_head.first].find_group(cluster_head.second) == cluster_head.second) {
                            result.erase(make_pair(cluster_head.first, combined_read_right[cluster_head.first]));
                            combined_read_right[cluster_head.first] = cluster_head.second;
                        } else {
                            result.erase(cluster_head);
                        }


                        tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].first = 
                            min_not_minus_one(best_left, best_right + loop_right + snarl_length);
                        tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].second = 
                            min_not_minus_one(best_right, loop_left == -1 ? -1 : best_left + loop_left + snarl_length);
                    }
                } else {
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first = 
                            min_not_minus_one(dist_left, dist_right + loop_right + snarl_length);
                }
            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "Found clusters on simple snarl " << snarl_index.id_in_parent << "-" << snarl_index.end_id << ": " << endl;
        for ( pair<size_t, size_t> cluster_head : result) {
            cerr << "\t" << tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).pos 
                << ": dist left: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first 
                << ", dist right: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second << endl;
        }
#endif

        return result;

    };
    
    ostream& operator<<(ostream& out, const SnarlSeedClusterer::Seed& seed) {
        return out << "Seed {"
            << "make_pos_t(" << id(seed.pos) << ", " << is_rev(seed.pos) << ", " << offset(seed.pos) << ")" << ", " 
            << seed.source << ", "
            << seed.is_top_level_node << ", "
            << seed.component << ", "
            << seed.offset << ", "
            << seed.is_top_level_snarl << ", "
            << seed.snarl_rank << ", "
            << seed.start_length << ", "
            << seed.end_length << ", "
            << seed.node_length << ", "
            << seed.rev_in_chain << "}"; 
    }
    
}
