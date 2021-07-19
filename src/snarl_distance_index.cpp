//#define debug_distance_indexing
//#define debug_snarl_traversal
//#define debug_distances

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {


make_distance_index(bdsg::SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit) {
    distance_index->set_snarl_size_limit(size_limit);

    //Build the temporary distance index from the graph
    TemporaryDistanceIndex temp_index(graph, snarl_finder, size_limit);

    //And fill in the permanent distance index
    vector<const TemporaryDistanceIndex*> indexes;
    indexes.emplace_back(&temp_index);
    get_snarl_tree_records(indexes, graph);
}

/*Fill in the snarl index.
 * The index will already know its boundaries and everything knows their relationships in the
 * snarl tree. This needs to fill in the distances and the ranks of children in the snarl
 * The rank of a child is arbitrary, except that the start node will always be 0 and the end node
 * will always be the node count+1 (since node count doesn't count the boundary nodes)
 */
void SnarlDistanceIndex::TemporaryDistanceIndex::populate_snarl_index(
                TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
                pair<temp_record_t, size_t> snarl_index, size_t size_limit, const HandleGraph* graph) {
#ifdef debug_distance_indexing
    cerr << "Getting the distances for snarl " << structure_start_end_as_string(snarl_index) << endl;
    assert(snarl_index.first == TEMP_SNARL);
#endif




    /*Helper function to find the ancestor of a node that is a child of this snarl */
    auto get_ancestor_of_node = [&](pair<temp_record_t, size_t> curr_index) {

        //This is a child that isn't a node, so it must be a chain
        if (curr_index.second == temp_snarl_record.start_node_id || 
            curr_index.second == temp_snarl_record.end_node_id) {
            return curr_index;
        }

        //Otherwise, walk up until we hit the current snarl
        pair<temp_record_t, size_t> parent_index = temp_node_records.at(curr_index.second-min_node_id).parent;
        while (parent_index != snarl_index) {
            curr_index=parent_index;
            parent_index = parent_index.first == TEMP_SNARL ? temp_snarl_records.at(parent_index.second).parent
                                                            : temp_chain_records.at(parent_index.second).parent;
#ifdef debug_distance_indexing
            assert(parent_index.first != TEMP_ROOT); 
#endif
        }
        
        return curr_index;
    };


    /*Now go through each of the children and add distances from that child to everything reachable from it
     * Start a dijkstra traversal from each node side in the snarl and record all distances
     */

    //Add the start and end nodes to the list of children so that we include them in the traversal 
    //TODO: Copying the list
    vector<pair<temp_record_t, size_t>> all_children = temp_snarl_record.children;
    if (!temp_snarl_record.is_root_snarl) {


        all_children.emplace_back(TEMP_NODE, temp_snarl_record.start_node_id);
        all_children.emplace_back(TEMP_NODE, temp_snarl_record.end_node_id);
    }

    while (!all_children.empty()) {
        const pair<temp_record_t, size_t> start_index = std::move(all_children.back());
        all_children.pop_back();

        //Check if this node is a tip
        if ((start_index.first == TEMP_NODE 
             && start_index.second != temp_snarl_record.start_node_id 
             && start_index.second != temp_snarl_record.end_node_id) 
            || 
            (start_index.first == TEMP_CHAIN && temp_chain_records.at(start_index.second).is_trivial)) {
            id_t node_id = start_index.first == TEMP_NODE ? start_index.second : temp_chain_records.at(start_index.second).start_node_id;
            size_t rank = start_index.first == TEMP_NODE ? temp_node_records.at(start_index.second-min_node_id).rank_in_parent 
                                                          : temp_chain_records.at(start_index.second).rank_in_parent;
            
            bool has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, false), false, [&](const handle_t next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_node_records.at(node_id-min_node_id).is_tip = true;
                temp_snarl_record.tippy_child_ranks.insert(rank);
            }
            has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, true), false, [&](const handle_t next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_node_records.at(node_id-min_node_id).is_tip = true;
                temp_snarl_record.tippy_child_ranks.insert(rank);
            }
        }

        bool start_is_tip = start_index.first == TEMP_NODE 
                      ? temp_node_records.at(start_index.second-min_node_id).is_tip 
                      : temp_chain_records.at(start_index.second).is_tip;

        size_t start_rank = start_index.first == TEMP_NODE 
                ? temp_node_records.at(start_index.second-min_node_id).rank_in_parent
                : temp_chain_records.at(start_index.second).rank_in_parent;


        if (start_index.first == TEMP_NODE && start_index.second == temp_snarl_record.start_node_id) {
            start_rank = 0;
        } else if (start_index.first == TEMP_NODE && start_index.second == temp_snarl_record.end_node_id) {
            start_rank = 1;
        } //TODO:
          //else {
          //  assert(start_rank != 0 && start_rank != 1);
          //}

        if ( (temp_snarl_record.node_count < size_limit || size_limit == 0) && !start_is_tip &&
             !start_rank == 0 && ! start_rank == 1) {
            //If we don't care about internal distances, and we also are not at a boundary or tip
            continue;
        }

        //Start from either direction for all nodes, but only going in for start and end
        vector<bool> directions;
        if (start_index.first == TEMP_NODE && start_index.second == temp_snarl_record.start_node_id) {
            directions.emplace_back(temp_snarl_record.start_node_rev);
        } else if (start_index.first == TEMP_NODE && start_index.second == temp_snarl_record.end_node_id){
            directions.emplace_back(!temp_snarl_record.end_node_rev);
        } else {
            directions.emplace_back(true);
            directions.emplace_back(false);
        }
        for (bool start_rev : directions) {
            //Start a dijkstra traversal from start_index going in the direction indicated by start_rev
            //Record the distances to each node (child of the snarl) found

#ifdef debug_distance_indexing
            cerr << "  Starting from child " << structure_start_end_as_string(start_index)
                 << " going " << (start_rev ? "rev" : "fd") << endl;
#endif

            //Define a NetgraphNode as the value for the priority queue:
            // <distance, <<type of node, index into temp_node/chain_records>, direction>
            using NetgraphNode = pair<int64_t, pair<pair<temp_record_t, size_t>, bool>>; 
            auto cmp = [] (const NetgraphNode a, const NetgraphNode b) {
                return a.first > b.first;
            };
            std::priority_queue<NetgraphNode, vector<NetgraphNode>, decltype(cmp)> queue(cmp);
            unordered_set<pair<pair<temp_record_t, size_t>, bool>> seen_nodes;
            queue.push(make_pair(0, make_pair(start_index, start_rev)));

            while (!queue.empty()) {

                int64_t current_distance = queue.top().first;
                pair<temp_record_t, size_t> current_index = queue.top().second.first;
                bool current_rev = queue.top().second.second;
                seen_nodes.emplace(queue.top().second);
                queue.pop();

                //The handle that we need to follow to get the next reachable nodes
                //If the current node is a node, then its just the node. Otherwise, it's the 
                //opposite side of the child chain
                handle_t current_end_handle = current_index.first == TEMP_NODE ? 
                        graph->get_handle(current_index.second, current_rev) :
                        (current_rev ? graph->get_handle(temp_chain_records[current_index.second].start_node_id, 
                                                        !temp_chain_records[current_index.second].start_node_rev) 
                                  : graph->get_handle(temp_chain_records[current_index.second].end_node_id, 
                                                      temp_chain_records[current_index.second].end_node_rev));
#ifdef debug_distance_indexing
                        cerr << "    at child " << structure_start_end_as_string(current_index) << " going "
                             << (current_rev ? "rev" : "fd") << " at actual node " << graph->get_id(current_end_handle) 
                             << (graph->get_is_reverse(current_end_handle) ? "rev" : "fd") << endl;
#endif
                graph->follow_edges(current_end_handle, false, [&](const handle_t next_handle) {
                    //At each of the nodes reachable from the current one, fill in the distance from the start
                    //node to the next node (current_distance). If this handle isn't leaving the snarl,
                    //add the next nodes along with the distance to the end of the next node
                    auto& node_record = temp_node_records.at(graph->get_id(next_handle)-min_node_id);
//TODO: The snarl decomposition should find tips now
//                    if (node_record.node_id == 0) {
//#ifdef debug_distance_indexing
//                        cerr << "Adding a tip " <<  graph->get_id(next_handle) << endl;
//#endif
//                        //If we haven't seen this node before, it means that it was a tip
//                        node_record.node_id = graph->get_id(next_handle);
//                        node_record.node_length = graph->get_length(next_handle);
//                        node_record.rank_in_parent = temp_snarl_record.node_count+2;
//                        node_record.reversed_in_parent = false;
//                        node_record.parent = snarl_index; 
//                        node_record.is_tip = true;
//
//                        //also update the parent
//                        temp_snarl_record.node_count ++;
//                        temp_snarl_record.is_trivial = false;
//                        temp_snarl_record.children.emplace_back(TEMP_NODE, graph->get_id(next_handle));
//                        temp_snarl_record.tippy_child_ranks.insert(node_record.rank_in_parent);
//
//                        //TODO: Is it bad to change the list as we're walking through it?
//                        all_children.emplace_back(TEMP_NODE, graph->get_id(next_handle)); 
//                    }

                    //The index of the snarl's child that next_handle represents
                    pair<temp_record_t, size_t> next_index = get_ancestor_of_node(make_pair(TEMP_NODE, graph->get_id(next_handle))); 

                    bool next_is_tip = start_index.first == TEMP_NODE 
                              ? temp_node_records.at(start_index.second-min_node_id).is_tip 
                              : temp_chain_records.at(start_index.second).is_tip;

                    //The rank and orientation of next in the snarl
                    size_t next_rank = next_index.first == TEMP_NODE 
                            ? node_record.rank_in_parent
                            : temp_chain_records[next_index.second].rank_in_parent;
                    if (next_index.first == TEMP_NODE && next_index.second == temp_snarl_record.start_node_id) {
                        next_rank = 0;
                    } else if (next_index.first == TEMP_NODE && next_index.second == temp_snarl_record.end_node_id) {
                        next_rank = 1;
                    } //TODO: This won't be true of root snarls 
                      //else {
                      //  assert(next_rank != 0 && next_rank != 1);
                      //}
                    bool next_rev = next_index.first == TEMP_NODE || temp_chain_records[next_index.second].is_trivial 
                            ? graph->get_is_reverse(next_handle) 
                            : graph->get_id(next_handle) == temp_chain_records[next_index.second].end_node_id;

                    if (size_limit != 0 &&
                        (temp_snarl_record.node_count < size_limit ||
                         (start_rank == 0 || start_rank == 1 || next_rank == 0 || next_rank == 1))) {
                        //If we are looking at all distances or we are looking at tips or boundaries

                        //Set the distance
                        pair<size_t, bool> start = !temp_snarl_record.is_root_snarl && (start_rank == 0 || start_rank == 1) 
                            ? make_pair(start_rank, false) : make_pair(start_rank, !start_rev);
                        pair<size_t, bool> next = !temp_snarl_record.is_root_snarl && (next_rank == 0 || next_rank == 1) 
                            ? make_pair(next_rank, false) : make_pair(next_rank, next_rev);
                        if (!temp_snarl_record.distances.count(make_pair(start, next)) ) {

                            temp_snarl_record.distances[make_pair(start, next)] = current_distance;
#ifdef debug_distance_indexing
                            cerr << "           Adding distance between ranks " << start.first << " " << start.second << " and " << next.first << " " << next.second << ": " << current_distance << endl;
#endif
                        }
                    }


                    if (seen_nodes.count(make_pair(next_index, next_rev)) == 0 &&
                        graph->get_id(next_handle) != temp_snarl_record.start_node_id &&
                        graph->get_id(next_handle) != temp_snarl_record.end_node_id) {
                        //If this isn't leaving the snarl, then add the next node to the queue, 
                        //along with the distance to traverse it
                        int64_t next_node_len = next_index.first == TEMP_NODE ? graph->get_length(next_handle) :
                                        temp_chain_records[next_index.second].min_length;
                        queue.push(make_pair(current_distance + next_node_len, 
                                       make_pair(next_index, next_rev)));
                    }
                    if (next_index.first == TEMP_CHAIN) {
                        int64_t loop_distance = next_rev ? temp_chain_records[next_index.second].backward_loops.back() 
                                                         : temp_chain_records[next_index.second].forward_loops.front();
                        if (loop_distance != std::numeric_limits<int64_t>::max() &&
                            seen_nodes.count(make_pair(next_index, !next_rev)) == 0 &&
                            graph->get_id(next_handle) != temp_snarl_record.start_node_id &&
                            graph->get_id(next_handle) != temp_snarl_record.end_node_id) {
                            //If the next node can loop back on itself, then add the next node in the opposite direction
                            int64_t next_node_len = loop_distance + 2 * graph->get_length(next_handle);
                            queue.push(make_pair(current_distance + next_node_len, 
                                           make_pair(next_index, !next_rev)));
                        }
                    }
#ifdef debug_distance_indexing
                    cerr << "        reached child " << structure_start_end_as_string(next_index) << "going " 
                         << (next_rev ? "rev" : "fd") << " with distance " << current_distance << " for ranks " << start_rank << " " << next_rank << endl;
#endif
                });
            }
        }
    }
}

void SnarlDistanceIndex::get_snarl_tree_records(const vector<const TemporaryDistanceIndex*>& temporary_indexes, const HandleGraph* graph) {

#ifdef debug_distance_indexing
    cerr << "Convert a temporary distance index into a permanent one" << endl;
#endif

    //TODO: Make sure not to include trivial chains
    //Convert temporary distance indexes into the final index stored as a single vector
    size_t total_index_size = 1;
    size_t total_component_count = 0;
    id_t min_node_id = 0;
    id_t max_node_id = 0;
    //The maximum value that gets stored in the index 
    //Used for setting the width of ints in the int vector (snarl_tree_records)
    //This only counts distance values, since the max offset value will be the length 
    //of the array, I'll check that right before setting the width
    size_t max_value = 0; 

    /*Go through each of the indexes to count how many nodes, components, etc */
    for (const TemporaryDistanceIndex* temp_index : temporary_indexes) {
        total_index_size += temp_index->index_size;
        total_component_count += temp_index->root_structure_count;
        min_node_id = min_node_id == 0 ? temp_index->min_node_id 
                                       : std::min(min_node_id, temp_index->min_node_id);
        max_node_id = std::max(max_node_id, temp_index->max_node_id);
    }

#ifdef debug_distance_indexing
    cerr << "Converting " << temporary_indexes.size() << " temporary indexes with "
         << total_component_count << " connected components from node " 
         << min_node_id << " to " << max_node_id << endl;
    cerr << " Adding root record" << endl;
#endif

    //TODO: Count everything properly
    //snarl_tree_records->reserve(total_index_size);

    /*Allocate memory for the root and the nodes */
    RootRecordConstructor root_record(0, total_component_count, max_node_id-min_node_id+1, min_node_id, &snarl_tree_records);
#ifdef debug_distance_indexing
    cerr << "  Root record had length " << snarl_tree_records->size() << endl;
#endif

    /*Now go through each of the chain/snarl indexes and copy them into snarl_tree_records
     * Walk down the snarl tree and fill in children
     */
    //TODO: For now I'm assuming that I'm including distances
    // maps <index into temporary_indexes, <record type, index into chain/snarl/node records>> to new offset
    std::unordered_map<pair<size_t, pair<temp_record_t, size_t>>, size_t> record_to_offset;
    //Set the root index
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        //Any root will point to the same root
        record_to_offset.emplace(make_pair(temp_index_i,make_pair(TEMP_ROOT, 0)), 0);
    }
    //Go through each separate temporary index, corresponding to separate connected components
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        const TemporaryDistanceIndex* temp_index = temporary_indexes[temp_index_i];
        //Get a stack of temporary snarl tree records to be added to the index
        //Initially, it contains only the root components
        //This reverses the order of the connected components but I don't think that matters
        //TODO: this is copying the components but it shouldn't be too big so I think it's fine
        vector<pair<temp_record_t, size_t>> temp_record_stack = temp_index->components;

        while (!temp_record_stack.empty()) {
            pair<temp_record_t, size_t> current_record_index = temp_record_stack.back();
            temp_record_stack.pop_back();

#ifdef debug_distance_indexing
            cerr << "Translating " << temp_index->structure_start_end_as_string(current_record_index) << endl;
#endif

            if (current_record_index.first == TEMP_CHAIN) {
                /*Add a new chain to the index. Each of the chain's child snarls and nodes will also 
                 * be added here
                 */
                const TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = 
                        temp_index->temp_chain_records[current_record_index.second];
                if (!temp_chain_record.is_trivial) {
                    //If this chain contains at least two nodes
#ifdef debug_distance_indexing
                    cerr << "  Adding this chain at offset " << snarl_tree_records->size() << endl;
                    cerr << "            with indices " << current_record_index.first << " " << current_record_index.second << endl;
#endif 
                    record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records->size());

                    ChainRecordConstructor chain_record_constructor;

                    if (temp_chain_record.chain_components.back() == 0 || snarl_size_limit == 0) {
                        record_t record_type = snarl_size_limit == 0 ? CHAIN : DISTANCED_CHAIN;
                        chain_record_constructor = ChainRecordConstructor(snarl_tree_records->size(), record_type, 
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                        chain_record_constructor.set_start_end_connected();
                    } else {
                        chain_record_constructor = ChainRecordConstructor(snarl_tree_records->size(), MULTICOMPONENT_CHAIN, 
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                    }
                    chain_record_constructor.set_parent_record_offset(
                            record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);//TODO: Get the actual parent

                    chain_record_constructor.set_min_length(temp_chain_record.min_length);
                    chain_record_constructor.set_max_length(temp_chain_record.max_length);
                    chain_record_constructor.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    chain_record_constructor.set_start_node(temp_chain_record.start_node_id, temp_chain_record.start_node_rev);
                    chain_record_constructor.set_end_node(temp_chain_record.end_node_id, temp_chain_record.end_node_rev);
                    max_value = std::max(max_value, (size_t) temp_chain_record.max_length);


                    size_t chain_node_i = 0; //How far along the chain are we?
                    bool prev_node = false;//Was the previous thing in the chain a node?
                    pair<size_t, bool> last_child_offset;

                    for (size_t child_record_index_i = 0 ; child_record_index_i < temp_chain_record.children.size() ; child_record_index_i++) {
                        const pair<temp_record_t, size_t>& child_record_index = temp_chain_record.children[child_record_index_i];
                        //Go through each node and snarl in the chain and add them to the index
#ifdef debug_distance_indexing
                        cerr << "  Adding child of the chain: " << temp_index->structure_start_end_as_string(child_record_index) << endl;
#endif

                        if (child_record_index.first == TEMP_NODE) {
                            //Add a node to the chain
                            if (prev_node) {
                                //If the last thing we saw was a node, then this is the end of a trivial snarl 
                                chain_record_constructor.add_trivial_snarl();
#ifdef debug_distance_indexing
                                cerr << "    Also adding a trivial snarl before this node" << endl;
#endif
                            }
                            if (chain_node_i != 0 && child_record_index == temp_chain_record.children.front()) {
                                //If this is the last node in the chain, and it is the same as the first node -
                                // it is a looping chain and we set this and don't re-record the node
                                // TODO: I'm using externally_start_start_connected here to indicate that 
                                // it's sharing a start and end node, but chains might actually be allowed to
                                // be start-end connected in which case I need a new flag
                                //chain_record_constructor.set_externally_start_end_connected();
#ifdef debug_distance_indexing
                            cerr << "    This is a looping chain"  << endl;
#endif

                            } else {
                                const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                                    temp_index->temp_node_records[child_record_index.second-min_node_id];
                                //Fill in this node in the index 
                                record_t record_type = snarl_size_limit == 0 ? NODE : DISTANCED_NODE;
                                NodeRecordConstructor node_record_constructor(
                                    get_offset_from_node_id(temp_node_record.node_id), record_type, &snarl_tree_records);
                                node_record_constructor.set_node_length(temp_node_record.node_length);
                                node_record_constructor.set_parent_record_offset(chain_record_constructor.get_offset());
                                node_record_constructor.set_is_reversed_in_parent(temp_node_record.reversed_in_parent);
                                max_value = std::max(max_value, temp_node_record.node_length);

                                //TODO: This is not really used
                                //The "rank" of the node actually points to the node in the chain, so it is the
                                //current size of the records (before adding the node to the chain)
                                node_record_constructor.set_rank_in_parent(snarl_tree_records->size());

                                //Add the node to the chain
                                if (chain_record_constructor.ChainRecord::get_record_type() == MULTICOMPONENT_CHAIN) {
                                    chain_record_constructor.add_node(temp_node_record.node_id, temp_chain_record.prefix_sum[chain_node_i],
                                                                  temp_chain_record.forward_loops[chain_node_i],
                                                                  temp_chain_record.backward_loops[chain_node_i],
                                                                  temp_chain_record.chain_components[chain_node_i]);
                                } else {
                                    chain_record_constructor.add_node(temp_node_record.node_id, temp_chain_record.prefix_sum[chain_node_i],
                                                                  temp_chain_record.forward_loops[chain_node_i],
                                                                  temp_chain_record.backward_loops[chain_node_i]);
                                }
                                last_child_offset = make_pair(node_record_constructor.NodeRecord::record_offset, false);


#ifdef debug_distance_indexing
                            cerr << "    The node record is at offset " << node_record_constructor.NodeRecord::record_offset << endl;
#endif
                            }

                            chain_node_i++;
                            prev_node = true;


                        } else {
                            //Add a snarl to the chain
                            assert(child_record_index.first == TEMP_SNARL);
                            //Get the temporary snarl record
                            const TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = 
                                 temp_index->temp_snarl_records[child_record_index.second];
                            if (!temp_snarl_record.is_trivial) {
                                //If this is an actual snarl that we need to make

                                //Record how to find the new snarl record
                                record_to_offset.emplace(make_pair(temp_index_i, child_record_index), snarl_tree_records->size()+1);

                                //Add the snarl to the chain, and get back the record to fill it in

                                record_t record_type = snarl_size_limit == 0 ? SNARL : 
                                    (temp_snarl_record.node_count < snarl_size_limit ? DISTANCED_SNARL : OVERSIZED_SNARL);
                                SnarlRecordConstructor snarl_record_constructor = 
                                    chain_record_constructor.add_snarl(temp_snarl_record.node_count, record_type);

                                //Fill in snarl info
                                snarl_record_constructor.set_parent_record_offset(chain_record_constructor.get_offset());
                                snarl_record_constructor.set_start_node(temp_snarl_record.start_node_id,
                                                                         temp_snarl_record.start_node_rev);
                                snarl_record_constructor.set_end_node(temp_snarl_record.end_node_id,
                                                                         temp_snarl_record.end_node_rev);
                                snarl_record_constructor.set_min_length(temp_snarl_record.min_length);
                                snarl_record_constructor.set_max_length(temp_snarl_record.max_length);
                                max_value = std::max(max_value, (size_t) temp_snarl_record.max_length);

                                //Add distances and record connectivity
                                for (const auto& it : temp_snarl_record.distances) {
                                    const pair<pair<size_t, bool>, pair<size_t, bool>>& node_ranks = it.first;
                                    const int64_t distance = it.second;

                                    if (snarl_size_limit != 0 &&
                                        (temp_snarl_record.node_count < snarl_size_limit ||
                                         (node_ranks.first.first == 0 || node_ranks.first.first == 1 ||
                                          node_ranks.second.first == 0 || node_ranks.second.first == 1))) {
                                        //If we are keeping track of distances and either this is a small enough snarl,
                                        //or the snarl is too big but we are looking at the boundaries
                                        snarl_record_constructor.set_distance(node_ranks.first.first, node_ranks.first.second, 
                                            node_ranks.second.first, node_ranks.second.second, distance);
                                        max_value = std::max(max_value, (size_t)distance);
                                    }

                                    //Now set the connectivity of this snarl
                                    if (node_ranks.first.first == 0 && node_ranks.second.first == 0) {
                                        snarl_record_constructor.set_start_start_connected();
                                    } else if ((node_ranks.first.first == 0 && node_ranks.second.first == 1) ||
                                               (node_ranks.first.first == 1 && node_ranks.second.first == 0)) {
                                        snarl_record_constructor.set_start_end_connected();
                                    } else if (node_ranks.first.first == 1 && node_ranks.second.first == 1) {
                                        snarl_record_constructor.set_end_end_connected();
                                    } else if ((node_ranks.first.first == 0 || node_ranks.second.first == 0) &&
                                               (temp_snarl_record.tippy_child_ranks.count(node_ranks.first.first) 
                                                || temp_snarl_record.tippy_child_ranks.count(node_ranks.second.first))) {
                                        snarl_record_constructor.set_start_tip_connected();
                                    } else if ((node_ranks.first.first == 1 || node_ranks.second.first == 1) &&
                                               (temp_snarl_record.tippy_child_ranks.count(node_ranks.first.first) 
                                                || temp_snarl_record.tippy_child_ranks.count(node_ranks.second.first))) {
                                        snarl_record_constructor.set_end_tip_connected();
                                    } else if (temp_snarl_record.tippy_child_ranks.count(node_ranks.first.first) 
                                                && temp_snarl_record.tippy_child_ranks.count(node_ranks.second.first)) {
                                        snarl_record_constructor.set_tip_tip_connected();
                                    }
                                }
#ifdef debug_distance_indexing
                            cerr << "    The snarl record is at offset " << snarl_record_constructor.SnarlRecord::record_offset << endl;
                            cerr << "    This child snarl has " << snarl_record_constructor.get_node_count() << " children: " << endl;
#endif
                                for (const pair<temp_record_t, size_t>& child : temp_snarl_record.children) {
                                    temp_record_stack.emplace_back(child);
#ifdef debug_distance_indexing
                                    cerr << "      " << temp_index->structure_start_end_as_string(child) << endl;
#endif
                                }


                                //Add connectivity in chain based on this snarl
                                //TODO: Tip-tip connected?
                                //TODO: What about if it's a multicomponent chain?
                                if (!snarl_record_constructor.get_is_reversed_in_parent()) {
                                    //If this snarl is oriented forward in the chain
                                    if (snarl_record_constructor.is_start_start_connected()) {
                                        chain_record_constructor.set_start_start_connected();
                                    } 
                                    if (snarl_record_constructor.is_end_end_connected()) {
                                        chain_record_constructor.set_end_end_connected();
                                    } 
                                    if (snarl_record_constructor.is_start_tip_connected()) {
                                        chain_record_constructor.set_start_tip_connected();
                                    }
                                    if (snarl_record_constructor.is_end_tip_connected()) {
                                        chain_record_constructor.set_end_tip_connected();
                                    }
                                } else {
                                    //If this snarl is oriented backwards in the chain
                                    if (snarl_record_constructor.is_start_start_connected()) {
                                        chain_record_constructor.set_end_end_connected();
                                    } 
                                    if (snarl_record_constructor.is_end_end_connected()) {
                                        chain_record_constructor.set_start_start_connected();
                                    } 
                                    if (snarl_record_constructor.is_start_tip_connected()) {
                                        chain_record_constructor.set_end_tip_connected();
                                    }
                                    if (snarl_record_constructor.is_end_tip_connected()) {
                                        chain_record_constructor.set_start_tip_connected();
                                    }
                                }
                                last_child_offset = make_pair(snarl_record_constructor.SnarlRecord::record_offset, true);
                            } else {
                                //Add a trivial snarl
#ifdef debug_distance_indexing
                                cerr << "    which is a trivial snarl" << endl;
#endif
                                chain_record_constructor.add_trivial_snarl();
                            }

                            prev_node = false;
                        }
                    }
                    //Does the chain loop and is the last node connected to the rest of the chain through the last snarl
                    bool last_node_connected = temp_chain_record.loopable && (temp_chain_record.start_node_id==temp_chain_record.end_node_id);
                    chain_record_constructor.set_last_child_offset(last_child_offset.first, last_child_offset.second, last_node_connected);
                    //Finish the chain by adding two 0's
                    chain_record_constructor.finish_chain();
                } else {
                    //If the chain is trivial, then only record the node
#ifdef debug_distance_indexing
                    cerr << "        this chain is actually just a node "
                         << temp_index->structure_start_end_as_string(temp_chain_record.children[0]) << endl; 
#endif
                    assert(temp_chain_record.children.size() == 1);
                    assert(temp_chain_record.children[0].first == TEMP_NODE);
                    const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                            temp_index->temp_node_records[temp_chain_record.children[0].second-min_node_id];

                    record_t record_type = snarl_size_limit == 0 ? NODE : DISTANCED_NODE;
                    NodeRecordConstructor node_record(
                        get_offset_from_node_id(temp_node_record.node_id), record_type, &snarl_tree_records);
                    node_record.set_node_length(temp_node_record.node_length);
                    node_record.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    node_record.set_is_reversed_in_parent(false);//TODO: Make sure that this is true
                    node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);
                    max_value = std::max(max_value, temp_node_record.node_length);

                    record_to_offset.emplace(make_pair(temp_index_i, current_record_index), node_record.NodeRecord::record_offset);

                }
            } else if (current_record_index.first == TEMP_SNARL) {
#ifdef debug_distance_indexing
                cerr << "        this is a root-level snarl " 
                     << temp_index->structure_start_end_as_string(current_record_index) << endl;
#endif
                //This is a root-level snarl
                record_t record_type = snarl_size_limit == 0 ? ROOT_SNARL : DISTANCED_ROOT_SNARL;

                const TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index->temp_snarl_records[current_record_index.second];
                record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records->size());

                SnarlRecordConstructor snarl_record_constructor (temp_snarl_record.node_count, &snarl_tree_records, record_type);

                //Fill in snarl info
                snarl_record_constructor.set_parent_record_offset(0);
                //Add distances and record connectivity
                for (const auto& it : temp_snarl_record.distances) {
                    const pair<pair<size_t, bool>, pair<size_t, bool>>& node_ranks = it.first;
                    const int64_t distance = it.second;
                    if (snarl_size_limit != 0 ) {
                        //TODO: I"m checking this but also automatically making a distanced snarl
                        //If we are keeping track of distances and either this is a small enough snarl,
                        //or the snarl is too big but we are looking at the boundaries
                        snarl_record_constructor.set_distance(node_ranks.first.first, node_ranks.first.second, 
                            node_ranks.second.first, node_ranks.second.second, distance);
                        max_value = std::max(max_value, (size_t)distance);
                    }
                }
#ifdef debug_distance_indexing
                cerr << "    The snarl record is at offset " << snarl_record_constructor.SnarlRecord::record_offset << endl;
                cerr << "    This child snarl has " << snarl_record_constructor.get_node_count() << " children: " << endl;
#endif
                for (const pair<temp_record_t, size_t>& child : temp_snarl_record.children) {
                        temp_record_stack.emplace_back(child);
                }
                
            } else {
                //TODO: This was a node that was a tip, so it wasn't put in a chain by the temporary index
                assert(current_record_index.first == TEMP_NODE);
#ifdef debug_distance_indexing
                cerr << "        this just a node that is a tip " 
                     << temp_index->structure_start_end_as_string(current_record_index) << endl;
#endif
                const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                        temp_index->temp_node_records[current_record_index.second-min_node_id];
                record_t record_type = snarl_size_limit == 0 ? NODE : DISTANCED_NODE;
                NodeRecordConstructor node_record(
                    get_offset_from_node_id(temp_node_record.node_id), record_type, &snarl_tree_records);
                node_record.set_node_length(temp_node_record.node_length);
                node_record.set_rank_in_parent(temp_node_record.rank_in_parent);
                node_record.set_is_reversed_in_parent(false);//TODO: Make sure that this is true
                node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_node_record.parent)]);
                max_value = std::max(max_value, temp_node_record.node_length);

                record_to_offset.emplace(make_pair(temp_index_i, current_record_index), node_record.NodeRecord::record_offset);
            }
#ifdef debug_distance_indexing
            cerr << "Finished translating " << temp_index->structure_start_end_as_string(current_record_index) << endl;
#endif
        }
#ifdef debug_distance_indexing
        cerr << "Adding roots" << endl;
#endif

        for (size_t component_num = 0 ; component_num < temp_index->components.size() ; component_num++){
            const pair<temp_record_t, size_t>& component_index = temp_index->components[component_num];
            //Let the root record know that it has another root
            root_record.add_component(component_num,record_to_offset[make_pair(temp_index_i,component_index)]); 

            SnarlTreeRecord record (record_to_offset[make_pair(temp_index_i, component_index)], 
                                    &snarl_tree_records);
            SnarlTreeRecordConstructor record_constructor(record_to_offset[make_pair(temp_index_i, component_index)], 
                                                              &snarl_tree_records);
            if (record.get_record_type() != ROOT_SNARL && record.get_record_type() != DISTANCED_ROOT_SNARL) {
                //If this isn't a root snarl
                handle_t start_out = graph->get_handle(record.get_start_id(), !record.get_start_orientation());
                handle_t end_out = graph->get_handle(record.get_end_id(), record.get_end_orientation());
                handle_t start_in = graph->get_handle(record.get_start_id(), record.get_start_orientation());
                handle_t end_in = graph->get_handle(record.get_end_id(), !record.get_end_orientation());


                graph->follow_edges(start_out, false, [&](const handle_t& h) {
                    if (h == start_in) {
                        record_constructor.set_externally_start_start_connected();
                    } else if (h == end_in) {
                        record_constructor.set_externally_start_end_connected();
                    }
                    return true;
                });
                graph->follow_edges(end_out, false, [&](const handle_t& h) {
                    if (h == end_in) {
                        record_constructor.set_externally_end_end_connected();
                    } else if (h == start_in) {
                        record_constructor.set_externally_start_end_connected();
                    }
                    return true;
                });
            }
            

#ifdef debug_distance_indexing
            cerr << temp_index->structure_start_end_as_string(component_index) << endl;
            //assert(record.get_parent_record_offset() == 0);
#endif
        }
    }

    

#ifdef debug_distance_indexing
    cerr << "Now filling in children of each snarl" << endl;
    cerr << "The index currently has size " << snarl_tree_records->size() << endl;
#endif

    /* Now go through everything again and give everything children */
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        const TemporaryDistanceIndex* temp_index = temporary_indexes[temp_index_i];
        for (size_t temp_snarl_i = 0 ; temp_snarl_i < temp_index->temp_snarl_records.size() ; temp_snarl_i ++) {
            //Get the temporary index for this snarl
            const TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index->temp_snarl_records[temp_snarl_i];
            if (!temp_snarl_record.is_trivial) {
                //And a constructor for the permanent record, which we've already created
                SnarlRecordConstructor snarl_record_constructor (&snarl_tree_records,
                        record_to_offset[make_pair(temp_index_i, make_pair(TEMP_SNARL, temp_snarl_i))]); 
                //Now add the children and tell the record where to find them
                snarl_record_constructor.set_child_record_pointer(snarl_tree_records->size());
                for (pair<temp_record_t, size_t> child : temp_snarl_record.children) {
                    snarl_record_constructor.add_child(record_to_offset[make_pair(temp_index_i, child)]);
#ifdef debug_distance_indexing
                cerr << "       child " << temp_index->structure_start_end_as_string(child) << endl;
                cerr << "        " << child.first << " " << child.second << endl;
                cerr << "        Add child " << net_handle_as_string(get_net_handle(record_to_offset[make_pair(temp_index_i, child)], START_END)) 
                     << "     at offset " << record_to_offset[make_pair(temp_index_i, child)] 
                     << "     to child list at offset " << snarl_tree_records->size() << endl;
#endif
                }
            }
        }
    }

    //TODO: Instead of repacking, could figure out the max size from the temporary index 
    //and set the width to start
    //The size of the array can be stored as an offset in the array
    max_value = std::max(max_value, snarl_tree_records->size());
    size_t new_width = log(max_value-1) + 1;
    new_width = std::max(new_width, (size_t)13); //13 is the width of the tags
    snarl_tree_records->repack(new_width, snarl_tree_records->size());
}
//TODO: Also need to go the other way, from final index to temporary one for merging


}
