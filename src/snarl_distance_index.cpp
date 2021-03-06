//#define debug_distance_indexing
//#define debug_snarl_traversal

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {


///////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor
SnarlDistanceIndex::SnarlDistanceIndex() {}
SnarlDistanceIndex::~SnarlDistanceIndex() {}

SnarlDistanceIndex::SnarlDistanceIndex(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder){

    //Build the temporary distance index from the graph
    TemporaryDistanceIndex temp_index(graph, snarl_finder);

    //And fill in the permanent distance index
    vector<const TemporaryDistanceIndex*> indexes;
    indexes.emplace_back(&temp_index);
    get_snarl_tree_records(indexes);
}

/*Temporary distance index for constructing the index from the graph
 */
SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryDistanceIndex(){}

SnarlDistanceIndex::TemporaryDistanceIndex::~TemporaryDistanceIndex(){}

SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryDistanceIndex(
    const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder) :
    min_node_id(graph->min_node_id()), max_node_id(graph->max_node_id()) {

#ifdef debug_distance_indexing
    cerr << "Creating new distance index for nodes between " << graph->min_node_id() << " and " << graph->max_node_id() << endl;

#endif
    //Construct the distance index using the snarl decomposition
    //traverse_decomposition will visit all structures (including trivial snarls), calling
    //each of the given functions for the start and ends of the snarls and chains
    temp_node_records.resize(max_node_id-min_node_id+1);


    
    //Stores unfinished records, as type of record and offset into appropriate vector
    //(temp_node/snarl/chain_records)
    vector<pair<temp_record_t, size_t>> stack;


    /*Go through the decomposition top down and record the connectivity of the snarls and chains
     * Distances will be added later*/

    snarl_finder->traverse_decomposition(
    [&](handle_t chain_start_handle) {
        /*This gets called when a new chain is found, starting at the start handle going into chain
         * For the first node in a chain, create a chain record and fill in the first node.
         * Also add the first node record
         */
#ifdef debug_distance_indexing
        cerr << "  Starting new chain at " << graph->get_id(chain_start_handle) << (graph->get_is_reverse(chain_start_handle) ? " reverse" : " forward") << endl;
        //We shouldn't have seen this node before
        //assert(temp_node_records[graph->get_id(chain_start_handle)-min_node_id].node_id == 0);
#endif

        //Fill in node in chain
        stack.emplace_back(TEMP_CHAIN, temp_chain_records.size());
        id_t node_id = graph->get_id(chain_start_handle);
        temp_chain_records.emplace_back();
        auto& temp_chain = temp_chain_records.back();
        temp_chain.start_node_id = node_id; 
        temp_chain.start_node_rev = graph->get_is_reverse(chain_start_handle);
        temp_chain.children.emplace_back(TEMP_NODE, node_id);

        //And the node record itself
        auto& temp_node = temp_node_records.at(node_id-min_node_id);
        temp_node.node_id = node_id;
        temp_node.node_length = graph->get_length(chain_start_handle);
        temp_node.rank_in_parent = 0;
        temp_node.reversed_in_parent = graph->get_is_reverse(chain_start_handle);
        temp_node.parent = stack.back(); //The parent is this chain

    },
    [&](handle_t chain_end_handle) {
        /*This gets called at the end of a chain, facing out
         * Record the chain's end node. The node record itself would have been added as part of the snarl
         * Also record the chain's parent here
         */

        //Done with this chain
        pair<temp_record_t, size_t> chain_index = stack.back();
        stack.pop_back();

        assert(chain_index.first == TEMP_CHAIN);
        TemporaryChainRecord& temp_chain_record = temp_chain_records.at(chain_index.second);
        id_t node_id = graph->get_id(chain_end_handle);

        //Fill in node in chain
        temp_chain_record.end_node_id = node_id;
        temp_chain_record.end_node_rev = graph->get_is_reverse(chain_end_handle);

        if (stack.empty()) {
            //If this was the last thing on the stack, then this was a root
            temp_chain_record.parent = make_pair(TEMP_ROOT, 0);
            root_structure_count += 1;
            components.emplace_back(chain_index);
        } else {
            //The last thing on the stack is the parent of this chain, which must be a snarl
            temp_chain_record.parent = stack.back();
            temp_snarl_records[temp_chain_record.parent.second].children.emplace_back(chain_index);
        }

        if (temp_chain_record.children.size() == 1 && temp_chain_record.start_node_id == temp_chain_record.end_node_id) {
            temp_chain_record.is_trivial = true;
        }


#ifdef debug_distance_indexing
        cerr << "  Ending new " << (temp_chain_record.is_trivial ? "trivial " : "") <<  "chain " << structure_start_end_as_string(chain_index)
             << endl << "    that is a child of " << structure_start_end_as_string(temp_chain_record.parent) << endl;
#endif
    },
    [&](handle_t snarl_start_handle) {
        /*This gets called at the beginning of a new snarl facing in
         * Create a new snarl record and fill in the start node.
         * The node record would have been created as part of the chain, or as the end node
         * of the previous snarl
         */

#ifdef debug_distance_indexing
        cerr << "  Starting new snarl at " << graph->get_id(snarl_start_handle) << (graph->get_is_reverse(snarl_start_handle) ? " reverse" : " forward") << endl;
#endif
        stack.emplace_back(TEMP_SNARL, temp_snarl_records.size());
        temp_snarl_records.emplace_back();
        temp_snarl_records.back().start_node_id = graph->get_id(snarl_start_handle);
        temp_snarl_records.back().start_node_rev = graph->get_is_reverse(snarl_start_handle);
        temp_snarl_records.back().start_node_length = graph->get_length(snarl_start_handle);
    },
    [&](handle_t snarl_end_handle){
        /*This gets called at the end of the snarl facing out
         * Fill in the end node of the snarl, its parent, and record the snarl as a child of its 
         * parent chain
         * Also create a node record
         */
        pair<temp_record_t, size_t> snarl_index = stack.back();
        stack.pop_back();
        assert(snarl_index.first == TEMP_SNARL);
        assert(stack.back().first == TEMP_CHAIN);
        TemporarySnarlRecord& temp_snarl_record = temp_snarl_records[snarl_index.second];
        id_t node_id = graph->get_id(snarl_end_handle);

        //Record the end node in the snarl
        temp_snarl_record.end_node_id = node_id;
        temp_snarl_record.end_node_rev = graph->get_is_reverse(snarl_end_handle);
        temp_snarl_record.end_node_length = graph->get_length(snarl_end_handle);
        temp_snarl_record.node_count = temp_snarl_record.children.size();
        if (temp_snarl_record.children.size() == 0) {
            temp_snarl_record.is_trivial = true;
        }
        //Record the snarl as a child of its chain
        if (stack.empty()) {
            assert(false);
            //TODO: The snarl should always be the child of a chain
            //If this was the last thing on the stack, then this was a root
            //TODO: I'm not sure if this would get put into a chain or not
            temp_snarl_record.parent = make_pair(TEMP_ROOT, 0);
            root_structure_count += 1;
            components.emplace_back(snarl_index);
        } else {
            //This is the child of a chain
            assert(stack.back().first == TEMP_CHAIN);
            temp_snarl_record.parent = stack.back();
            auto& temp_chain = temp_chain_records.at(stack.back().second);
            temp_chain.children.emplace_back(snarl_index);
            temp_chain.children.emplace_back(TEMP_NODE, node_id);

        }

        //Record the node itself. This gets done for the start of the chain, and ends of snarls
        TemporaryNodeRecord& temp_node_record = temp_node_records.at(node_id-min_node_id);
        temp_node_record.node_id = node_id;
        temp_node_record.node_length = graph->get_length(snarl_end_handle);
        temp_node_record.reversed_in_parent = graph->get_is_reverse(snarl_end_handle);
        temp_node_record.parent = stack.back();
        
        //TODO: This isn't actually counting everything
        index_size += SnarlRecord::record_size(DISTANCED_SNARL, temp_snarl_record.node_count);

#ifdef debug_distance_indexing
        cerr << "  Ending new snarl " << structure_start_end_as_string(snarl_index)
             << endl << "    that is a child of " << structure_start_end_as_string(temp_snarl_record.parent) << endl;
#endif
    });


    /*Now go through the decomposition again to fill in the distances
     * This traverses all chains in reverse order that we found them in, so bottom up
     * Each chain and snarl already knows its parents and children, except for single nodes
     * that are children of snarls. These nodes were not in chains will have their node 
     * records created here
     * TODO: I don't think the decomposition would have visited single nodes
     */
#ifdef debug_distance_indexing
    cerr << "Filling in the distances in snarls" << endl;
#endif
    for (int i = temp_chain_records.size()-1 ; i >= 0 ; i--) {

        TemporaryChainRecord& temp_chain_record = temp_chain_records[i];
#ifdef debug_distance_indexing
        cerr << "  At "  << (temp_chain_record.is_trivial ? " trivial " : "") << " chain " << structure_start_end_as_string(make_pair(TEMP_CHAIN, i)) << endl; 
#endif

        //Add the first values for the prefix sum and backwards loop vectors
        temp_chain_record.prefix_sum.emplace_back(0);
        temp_chain_record.backward_loops.emplace_back(std::numeric_limits<int64_t>::max());


        /*First, go through each of the snarls in the chain in the forward direction and
         * fill in the distances in the snarl. Also fill in the prefix sum and backwards
         * loop vectors here
         */
        for (const pair<temp_record_t, size_t>& chain_child_index : temp_chain_record.children){ 
            //Go through each of the children in the chain, skipping nodes 
            //The snarl may be trivial, in which case don't fill in the distances
#ifdef debug_distance_indexing
            cerr << "    Looking at child " << structure_start_end_as_string(chain_child_index) << endl; 
#endif

            if (chain_child_index.first == TEMP_SNARL){
                //This is where all the work gets done. Need to go through the snarl and add 
                //all distances, then add distances to the chain that this is in
                //The parent chain will be the last thing in the stack
                TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(chain_child_index.second);

//TODO: Some snarls may appear to be trivial but have tips, so we still need to look for extra nodes
//                if (temp_snarl_record.is_trivial) {
//                    //For a trivial snarl, don't bother filling in the distances but still add the prefix sum vector
//                    temp_chain_record.prefix_sum.emplace_back(temp_chain_record.prefix_sum.back() + 
//                                                              temp_snarl_record.start_node_length);
//                    //TODO: I think that adding to max() should still be infinity???
//                    temp_chain_record.backward_loops.emplace_back(temp_chain_record.backward_loops.back() + 
//                                                                  2 * temp_snarl_record.start_node_length);
//                } else {

                    //Fill in this snarl's distances
                    populate_snarl_index(temp_snarl_record, chain_child_index, graph);

                    int64_t length;
                    //TODO: Double check these orientations
                    if (temp_snarl_record.distances.count(make_pair(make_pair(0, false), make_pair(1, false)))){
                        length = temp_snarl_record.distances.at(make_pair(make_pair(0, false), make_pair(1, false)));
                    } else if (temp_snarl_record.distances.count(make_pair(make_pair(1, true), make_pair(0, false)))){
                        length = temp_snarl_record.distances.at(make_pair(make_pair(1, true), make_pair(0, false)));
                    } else {
                        length = std::numeric_limits<int64_t>::max();
                    }
                    temp_snarl_record.min_length = length;

                    //And get the distance values for the end node in the chain
                    temp_chain_record.prefix_sum.emplace_back(temp_chain_record.prefix_sum.back() + 
                        temp_snarl_record.min_length + temp_snarl_record.start_node_length);
                    temp_chain_record.backward_loops.emplace_back(std::min(temp_snarl_record.loop_end,
                        temp_chain_record.backward_loops.back() 
                        + 2 * (temp_snarl_record.start_node_length + temp_snarl_record.min_length)));
//                }
            }
        }

        /*Now that we've gone through all the snarls in the chain, fill in the forward loop vector 
         * by going through the chain in the backwards direction
         */
        temp_chain_record.forward_loops.resize(temp_chain_record.children.size(), 
                                               std::numeric_limits<int64_t>::max());
        for (int j = (int)temp_chain_record.children.size() - 2 ; j >= 0 ; j--) {
            // We start at the next to last record because we need to look at this record and the next one.
            auto& child = temp_chain_record.children.at(j);
            if (child.first == TEMP_SNARL){
                TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(child.second);
                if (temp_snarl_record.is_trivial) {
                    temp_chain_record.forward_loops.at(j) = temp_chain_record.forward_loops.at(j+1) + 
                                    2*temp_snarl_record.end_node_length;
                } else {
                    temp_chain_record.forward_loops.at(j) = std::min(temp_chain_record.forward_loops.at(j+1) + 
                                  2*temp_snarl_record.end_node_length, temp_snarl_record.loop_start);
                }
            }
        }
    }
}

/*Fill in the snarl index.
 * The index will already know its boundaries and everything knows their relationships in the
 * snarl tree. This needs to fill in the distances and the ranks of children in the snarl
 * The rank of a child is arbitrary, except that the start node will always be 0 and the end node
 * will always be the node count+1 (since node count doesn't count the boundary nodes)
 */
void SnarlDistanceIndex::TemporaryDistanceIndex::populate_snarl_index(
                TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, 
                pair<temp_record_t, size_t> snarl_index, const HandleGraph* graph) {
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

    /*The boundary handles for this snarl */
    handle_t snarl_start_in = graph->get_handle(temp_snarl_record.start_node_id, temp_snarl_record.start_node_rev);
    handle_t snarl_start_out = graph->get_handle(temp_snarl_record.start_node_id, !temp_snarl_record.start_node_rev);
    handle_t snarl_end_out = graph->get_handle(temp_snarl_record.end_node_id, temp_snarl_record.end_node_rev);
    handle_t snarl_end_in = graph->get_handle(temp_snarl_record.end_node_id, !temp_snarl_record.end_node_rev);



    /*Now go through each of the children and add distances from that child to everything reachable from it
     * Start a dijkstra traversal from each node side in the snarl and record all distances
     */

    //Add the start and end nodes to the list of children so that we include them in the traversal 
    //TODO: Copying the list
    vector<pair<temp_record_t, size_t>> all_children = temp_snarl_record.children;
    all_children.emplace_back(TEMP_NODE, temp_snarl_record.start_node_id);
    all_children.emplace_back(TEMP_NODE, temp_snarl_record.end_node_id);

    for (const pair<temp_record_t, size_t>& start_index : all_children) {

        size_t start_rank = start_index.first == TEMP_NODE 
                ? (start_index.second == temp_snarl_record.start_node_id ? 0 : 1)
                : temp_chain_records.at(start_index.second).rank_in_parent;
        //Start from either direction for all nodes, but only going in for start and end
        vector<bool> directions;
        if (start_index.second == temp_snarl_record.start_node_id) {
            directions.emplace_back(temp_snarl_record.start_node_rev);
        } else if (start_index.second == temp_snarl_record.end_node_id){
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
                    if (node_record.node_id == 0) {
#ifdef debug_distance_indexing
                        cerr << "Adding a tip " <<  graph->get_id(next_handle) << endl;
#endif
                        //If we haven't seen this node before, it means that it was a tip
                        node_record.node_id = graph->get_id(next_handle);
                        node_record.node_length = graph->get_length(next_handle);
                        node_record.rank_in_parent = temp_snarl_record.node_count+2;
                        node_record.reversed_in_parent = false;
                        node_record.parent = snarl_index; 

                        //also update the parent
                        temp_snarl_record.node_count ++;
                        temp_snarl_record.is_trivial = false;
                        temp_snarl_record.children.emplace_back(TEMP_NODE, graph->get_id(next_handle));
                    }

                    //The index of the snarl's child that next_handle represents
                    pair<temp_record_t, size_t> next_index = get_ancestor_of_node(make_pair(TEMP_NODE, graph->get_id(next_handle))); 

                    //The rank and orientation of next in the snarl
                    size_t next_rank = next_index.first == TEMP_NODE 
                            ? node_record.rank_in_parent
                            : temp_chain_records[next_index.second].rank_in_parent;
                    bool next_rev = next_index.first == TEMP_NODE ? graph->get_is_reverse(next_handle) :
                            graph->get_id(next_handle) == temp_chain_records[next_index.second].end_node_id;

                    //Set the distance
                    temp_snarl_record.distances[make_pair(make_pair(start_rank, !start_rev), make_pair(next_rank, next_rev))] = current_distance;


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
#ifdef debug_distance_indexing
                    cerr << "        reached child " << structure_start_end_as_string(next_index) << "going " 
                         << (next_rev ? "rev" : "fd") << " with distance " << current_distance << endl;
#endif
                });
            }
        }
    }
}

string SnarlDistanceIndex::TemporaryDistanceIndex::structure_start_end_as_string(pair<temp_record_t, size_t> index) const {
    if (index.first == TEMP_NODE) {
        assert(index.second == temp_node_records[index.second-min_node_id].node_id);
        return "node " + std::to_string(temp_node_records[index.second-min_node_id].node_id);
    } else if (index.first == TEMP_SNARL) {
        const TemporarySnarlRecord& temp_snarl_record =  temp_snarl_records[index.second];
        return "snarl " + std::to_string( temp_snarl_record.start_node_id) 
                + (temp_snarl_record.start_node_rev ? " rev" : " fd") 
                + " -> " + std::to_string( temp_snarl_record.end_node_id) 
                + (temp_snarl_record.end_node_rev ? " rev" : " fd");
    } else if (index.first == TEMP_CHAIN) {
        const TemporaryChainRecord& temp_chain_record = temp_chain_records[index.second];
        return "chain " + std::to_string( temp_chain_record.start_node_id) 
                + (temp_chain_record.start_node_rev ? " rev" : " fd") 
                + " -> "  + std::to_string( temp_chain_record.end_node_id) 
                + (temp_chain_record.end_node_rev ? " rev" : " fd");
    } else if (index.first == TEMP_ROOT) {
        return (string) "root";
    } else {
        return (string)"???" + std::to_string(index.first) + "???";
    }
}

vector<size_t> SnarlDistanceIndex::get_snarl_tree_records(const vector<const TemporaryDistanceIndex*>& temporary_indexes) {

#ifdef debug_distance_indexing
    cerr << "Convert a temporary distance index into a permanent one" << endl;
#endif

    //TODO: Make sure not to include trivial chains
    //Convert temporary distance indexes into the final index stored as a single vector
    size_t total_index_size = 1;
    size_t total_component_count = 0;
    id_t min_node_id = 0;
    id_t max_node_id = 0;

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
    //snarl_tree_records.reserve(total_index_size);

    /*Allocate memory for the root and the nodes */
    //TODO: Could also write directly into snarl_tree_records if I allocate the right amount of memory
    RootRecordConstructor root_record(0, total_component_count, max_node_id-min_node_id+1, min_node_id, &snarl_tree_records);
#ifdef debug_distance_indexing
    cerr << "  Root record had length " << snarl_tree_records.size() << endl;
#endif

    /*Now go through each of the chain/snarl indexes and copy them into snarl_tree_records
     * Walk down the snarl tree and fill in children
     */
    //TODO: For now I'm assuming that I'm including distances
    //TODO: What about connectivity?
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
                    cerr << "  Adding this chain at offset " << snarl_tree_records.size() << endl;
#endif 
                    record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records.size());

                    ChainRecordConstructor chain_record_constructor(snarl_tree_records.size(), DISTANCED_CHAIN, 
                                                                  temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                    chain_record_constructor.set_parent_record_offset(
                            record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);//TODO: Get the actual parent
                    chain_record_constructor.set_min_length(temp_chain_record.min_length);
                    chain_record_constructor.set_max_length(temp_chain_record.max_length);
                    chain_record_constructor.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    chain_record_constructor.set_start_node(temp_chain_record.start_node_id, temp_chain_record.start_node_rev);
                    chain_record_constructor.set_end_node(temp_chain_record.end_node_id, temp_chain_record.end_node_rev);
                    size_t chain_node_i = 0; //How far along the chain are we?
                    bool prev_node = false;//Was the previous thing in the chain a node?

                    for (const pair<temp_record_t, size_t>& child_record_index : temp_chain_record.children) {
                        //Go through each node and snarl in the chain and add them to the index
#ifdef debug_distance_indexing
                        cerr << "  Adding child of the chain: " << temp_index->structure_start_end_as_string(child_record_index) << endl;
#endif

                        if (child_record_index.first == TEMP_NODE) {
                            //Add a node to the chain
                            if (prev_node) {
                                //TODO: I think trivial snarls would actually get their own temp record
                                //If the last thing we saw was a node, then this is the end of a trivial snarl 
                                chain_record_constructor.add_trivial_snarl();
#ifdef debug_distance_indexing
                                cerr << "    Also adding a trivial snarl before this node" << endl;
#endif
                            }
                            const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                                temp_index->temp_node_records[child_record_index.second-min_node_id];
                            //Fill in this node in the index 
                            NodeRecordConstructor node_record_constructor(
                                get_offset_from_node_id(temp_node_record.node_id), DISTANCED_NODE, &snarl_tree_records);
                            node_record_constructor.set_node_length(temp_node_record.node_length);
                            node_record_constructor.set_parent_record_offset(chain_record_constructor.get_offset());
                            node_record_constructor.set_is_rev_in_parent(temp_node_record.reversed_in_parent);

                            //The "rank" of the node actually points to the node in the chain, so it is the
                            //current size of the records (before adding the node to the chain)
                            node_record_constructor.set_rank_in_parent(snarl_tree_records.size());

                            //Add the node to the chain
                            chain_record_constructor.add_node(temp_node_record.node_id, temp_chain_record.prefix_sum[chain_node_i],
                                                              temp_chain_record.forward_loops[chain_node_i],
                                                              temp_chain_record.backward_loops[chain_node_i]);

                            chain_node_i++;
                            prev_node = true;

#ifdef debug_distance_indexing
                            cerr << "    The node record is at offset " << node_record_constructor.NodeRecord::record_offset << endl;
#endif

                        } else {
                            //Add a snarl to the chain
                            assert(child_record_index.first == TEMP_SNARL);
                            //Get the temporary snarl record
                            const TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = 
                                 temp_index->temp_snarl_records[child_record_index.second];
                            if (!temp_snarl_record.is_trivial) {
                                //If this is an actual snarl that we need to make
                                record_to_offset.emplace(make_pair(temp_index_i, child_record_index), snarl_tree_records.size()+1);
                                //Add the snarl to the chain, and get back the record to fill it in
                                SnarlRecordConstructor snarl_record_constructor = 
                                    chain_record_constructor.add_snarl(temp_snarl_record.node_count, DISTANCED_SNARL);
                                snarl_record_constructor.set_parent_record_offset(chain_record_constructor.get_offset());
                                snarl_record_constructor.set_start_node(temp_snarl_record.start_node_id,
                                                                         temp_snarl_record.start_node_rev);
                                snarl_record_constructor.set_end_node(temp_snarl_record.end_node_id,
                                                                         temp_snarl_record.end_node_rev);
                                snarl_record_constructor.set_min_length(temp_snarl_record.min_length);
                                snarl_record_constructor.set_max_length(temp_snarl_record.max_length);
                                //TODO: Add distance
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
                    //Finish the chain by adding two 0's
                    snarl_tree_records.emplace_back(0);
                    snarl_tree_records.emplace_back(0);
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
                    NodeRecordConstructor node_record(
                        get_offset_from_node_id(temp_node_record.node_id), DISTANCED_NODE, &snarl_tree_records);
                    node_record.set_node_length(temp_node_record.node_length);
                    node_record.set_rank_in_parent(temp_node_record.rank_in_parent);
                    node_record.set_is_rev_in_parent(false);//TODO: Make sure that this is true
                    node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);

                    record_to_offset.emplace(make_pair(temp_index_i, current_record_index), node_record.NodeRecord::record_offset);
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
                    NodeRecordConstructor node_record(
                        get_offset_from_node_id(temp_node_record.node_id), DISTANCED_NODE, &snarl_tree_records);
                    node_record.set_node_length(temp_node_record.node_length);
                    node_record.set_rank_in_parent(temp_node_record.rank_in_parent);
                    node_record.set_is_rev_in_parent(false);//TODO: Make sure that this is true
                    node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_node_record.parent)]);

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
            root_record.add_component(component_num,record_to_offset[make_pair(temp_index_i,component_index)]); 
#ifdef debug_distance_indexing
            cerr << temp_index->structure_start_end_as_string(component_index) << endl;
#endif
        }
    }

    

#ifdef debug_distance_indexing
    cerr << "Now filling in children of each snarl" << endl;
    cerr << "The index currently has size " << snarl_tree_records.size() << endl;
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
                snarl_record_constructor.set_child_record_pointer(snarl_tree_records.size());
                for (pair<temp_record_t, size_t> child : temp_snarl_record.children) {
                    snarl_record_constructor.add_child(record_to_offset[make_pair(temp_index_i, child)]);
#ifdef debug_distance_indexing
                cerr << "    Add child " << net_handle_as_string(get_net_handle(record_to_offset[make_pair(temp_index_i, child)], START_END, CHAIN_HANDLE)) 
                     << " at offset " << record_to_offset[make_pair(temp_index_i, child)] 
                     << " to child list at offset " << snarl_tree_records.size() << endl;
#endif
            }
            }
        }
    }
    return snarl_tree_records;
}
//TODO: Also need to go the other way, from final index to temporary one for merging

///////////////////////////////////////////////////////////////////////////////////////////////////
//Implement the SnarlDecomposition's functions for moving around the snarl tree
//


net_handle_t SnarlDistanceIndex::get_root() const {
    // The root is the first thing in the index, the traversal is tip to tip
    return get_net_handle(0, START_END, ROOT_HANDLE);
}

bool SnarlDistanceIndex::is_root(const net_handle_t& net) const {
    return get_handle_type(net) == ROOT_HANDLE && get_record_offset(net) == 0;
}

bool SnarlDistanceIndex::is_snarl(const net_handle_t& net) const {
    return get_handle_type(net) == SNARL_HANDLE 
            && SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == SNARL_HANDLE;
}

bool SnarlDistanceIndex::is_chain(const net_handle_t& net) const {
    return get_handle_type(net) == CHAIN_HANDLE
            && (SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == CHAIN_HANDLE
            || SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == NODE_HANDLE);
}

bool SnarlDistanceIndex::is_trivial_chain(const net_handle_t& net) const {
    return get_handle_type(net) == CHAIN_HANDLE
            && SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == NODE_HANDLE;
}

bool SnarlDistanceIndex::is_node(const net_handle_t& net) const {
    return get_handle_type(net) == NODE_HANDLE 
            && SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == NODE_HANDLE;
}
bool SnarlDistanceIndex::is_sentinel(const net_handle_t& net) const {
    return get_handle_type(net) == SENTINEL_HANDLE
            && SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == SNARL_HANDLE;
}

net_handle_t SnarlDistanceIndex::get_net(const handle_t& handle, const handlegraph::HandleGraph* graph) const{
    return get_net_handle(get_offset_from_node_id(graph->get_id(handle)), 
                          graph->get_is_reverse(handle) ? END_START : START_END, 
                          NODE_HANDLE);
}
handle_t SnarlDistanceIndex::get_handle(const net_handle_t& net, const handlegraph::HandleGraph* graph) const{
    //TODO: Maybe also want to be able to get the graph handle of a sentinel
    if (get_handle_type(net) == SENTINEL_HANDLE) {
        SnarlRecord snarl_record(net, &snarl_tree_records);
        if (starts_at(net) == START) {
            return graph->get_handle(snarl_record.get_start_id(), 
                       ends_at(net) == START ? !snarl_record.get_start_orientation()   //Going out
                                             : snarl_record.get_start_orientation());  //Going in
        } else {
            assert (starts_at(net) == END);
            return graph->get_handle(snarl_record.get_end_id(), 
                       ends_at(net) == END ? snarl_record.get_end_orientation()   //Going out
                                           : !snarl_record.get_end_orientation());  //Going in
        }
    } else if (get_handle_type(net) == NODE_HANDLE || get_handle_type(net) == CHAIN_HANDLE ) {
        NodeRecord node_record(net, &snarl_tree_records);
        return graph->get_handle(get_node_id_from_offset(get_record_offset(net)), 
                                 get_connectivity(net) == START_END ? false : true);
    } else {
        throw runtime_error("error: trying to get a handle from a snarl, chain, or root");
    }
}

net_handle_t SnarlDistanceIndex::get_parent(const net_handle_t& child) const {

    //If the child is the sentinel of a snarl, just return the snarl
    if (get_handle_type(child) == SENTINEL_HANDLE) {
        return get_net_handle(get_record_offset(child), START_END, SNARL_HANDLE); 
    } else if (get_handle_type(child) == ROOT_HANDLE) {
        throw runtime_error("error: trying to find the parent of the root");
    } 

    //Otherwise, we need to move up one level in the snarl tree

    //Get the pointer to the parent, and keep the connectivity of the current handle
    size_t parent_pointer = SnarlTreeRecord(child, &snarl_tree_records).get_parent_record_offset();
    connectivity_t child_connectivity = get_connectivity(child);

    //TODO: I"m going into the parent record here, which could be avoided if things knew what their parents were, but I think if you're doing this you'd later go into the parent anyway so it's probably fine
    net_handle_record_t parent_type = SnarlTreeRecord(parent_pointer, &snarl_tree_records).get_record_handle_type();
    connectivity_t parent_connectivity = START_END;
    if ((child_connectivity == START_END || child_connectivity == END_START) 
        && (parent_type == CHAIN_HANDLE)) {
        //TODO: This also needs to take into account the orientation of the child, which I might be able to get around?
        parent_connectivity = child_connectivity;
    }
    if (get_handle_type(child) == NODE_HANDLE && parent_type != CHAIN_HANDLE) {
        //If this is a node and it's parent is not a chain, we want to pretend that its 
        //parent is a chain
        return get_net_handle(get_record_offset(child), parent_connectivity, CHAIN_HANDLE);
    }

    return get_net_handle(parent_pointer, parent_connectivity);
}

net_handle_t SnarlDistanceIndex::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    if (get_handle_type(snarl) == CHAIN_HANDLE) {
        id_t id = get_end ? SnarlTreeRecord(snarl, &snarl_tree_records).get_end_id() 
                          : SnarlTreeRecord(snarl, &snarl_tree_records).get_start_id();
        bool rev_in_parent = NodeRecord(get_offset_from_node_id(id), &snarl_tree_records).get_is_rev_in_parent();
        if (get_end) {
            rev_in_parent = !rev_in_parent;
        }
        if (!face_in){
            rev_in_parent = !rev_in_parent;
        }
        connectivity_t connectivity = rev_in_parent ? END_START : START_END;
        return get_net_handle(get_offset_from_node_id(id), connectivity,  NODE_HANDLE);
    } else {
        assert(get_handle_type(snarl) == SNARL_HANDLE);
        endpoint_t start = get_end ? END : START;
        endpoint_t end = face_in ? (start == END ? START : END) : start;
        return get_net_handle(get_record_offset(snarl), endpoints_to_connectivity(start, end), SENTINEL_HANDLE);
    }
}

net_handle_t SnarlDistanceIndex::flip(const net_handle_t& net) const {
    connectivity_t old_connectivity = get_connectivity(net);
    connectivity_t new_connectivity =  endpoints_to_connectivity(get_end_endpoint(old_connectivity), 
                                                                get_start_endpoint(old_connectivity));
    return get_net_handle(get_record_offset(net), new_connectivity, get_handle_type(net));
}

net_handle_t SnarlDistanceIndex::canonical(const net_handle_t& net) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    connectivity_t connectivity;
    if (record.is_start_end_connected()) {
        connectivity = START_END;
    } else if (record.is_start_tip_connected()) {
        connectivity = START_TIP;
    } else if (record.is_end_tip_connected()) {
        connectivity = END_TIP;
    } else if (record.is_start_start_connected()) {
        connectivity = START_START;
    } else if (record.is_end_end_connected()) {
        connectivity = END_END;
    } else if (record.is_tip_tip_connected()) {
        connectivity = TIP_TIP;
    } else {
        connectivity = START_END; //TODO: put this back throw runtime_error("error: This node has no connectivity");
    }
    return get_net_handle(get_record_offset(net), connectivity);
}

SnarlDecomposition::endpoint_t SnarlDistanceIndex::starts_at(const net_handle_t& traversal) const {
    return get_start_endpoint(get_connectivity(traversal));

}
SnarlDecomposition::endpoint_t SnarlDistanceIndex::ends_at(const net_handle_t& traversal) const {
    return get_end_endpoint( get_connectivity(traversal));
}

//TODO: I'm also allowing this for the root
bool SnarlDistanceIndex::for_each_child_impl(const net_handle_t& traversal, const std::function<bool(const net_handle_t&)>& iteratee) const {
#ifdef debug_snarl_traversal
    cerr << "Go through children of " << net_handle_as_string(traversal) << endl;
#endif
    //What is this according to the snarl tree
    net_handle_record_t record_type = SnarlTreeRecord(traversal, &snarl_tree_records).get_record_handle_type();
    //What is this according to the handle 
    //(could be a trivial chain but actually a node according to the snarl tree)
    net_handle_record_t handle_type = get_handle_type(traversal);
    if (record_type == SNARL_HANDLE) {
        SnarlRecord snarl_record(traversal, &snarl_tree_records);
        return snarl_record.for_each_child(iteratee);
    } else if (record_type == CHAIN_HANDLE) {
        ChainRecord chain_record(traversal, &snarl_tree_records);
        return chain_record.for_each_child(iteratee);
    } else if (record_type == ROOT_HANDLE) {
        RootRecord root_record(traversal, &snarl_tree_records);
        return root_record.for_each_child(iteratee);
    } else if (record_type == NODE_HANDLE && handle_type == CHAIN_HANDLE) {
        //This is actually a node but we're pretending it's a chain
#ifdef debug_snarl_traversal
        cerr << "     which is actually a node pretending to be a chain" << endl;
#endif
        return iteratee(get_net_handle(get_record_offset(traversal), get_connectivity(traversal), NODE_HANDLE));
    } else {
        throw runtime_error("error: Looking for children of a node or sentinel");
    }
   
}

bool SnarlDistanceIndex::for_each_traversal_impl(const net_handle_t& item, const std::function<bool(const net_handle_t&)>& iteratee) const {
    if (get_handle_type(item) == SENTINEL_HANDLE) {
        //TODO: I'm not sure what to do here?
        if (!iteratee(get_net_handle(get_record_offset(item), START_END, get_handle_type(item)))) {
            return false;
        }
        if (!iteratee(get_net_handle(get_record_offset(item), END_START, get_handle_type(item)))) {
            return false;
        }
    }
    SnarlTreeRecord record(item, &snarl_tree_records);
    for ( size_t type = 1 ; type <= 9 ; type ++ ){
        connectivity_t connectivity = static_cast<connectivity_t>(type);
        if (record.has_connectivity(connectivity)) {
            if (!iteratee(get_net_handle(get_record_offset(item), connectivity, get_handle_type(item)))) {
                return false;
            }
        }
    }
    return true;
}

bool SnarlDistanceIndex::follow_net_edges_impl(const net_handle_t& here, const handlegraph::HandleGraph* graph, bool go_left, const std::function<bool(const net_handle_t&)>& iteratee) const {
#ifdef debug_snarl_traversal
    cerr << "following edges from " << net_handle_as_string(here) << " going " << (go_left ? "rev" : "fd") << endl;
    cerr << "        that is a child of " << net_handle_as_string(get_parent(here)) << endl;
#endif

    SnarlTreeRecord this_record(here, &snarl_tree_records);
    SnarlTreeRecord parent_record (get_parent(here), &snarl_tree_records);

    if (parent_record.get_record_handle_type() == ROOT_HANDLE) {
        //TODO: I'm not sure what to do in this case
        throw runtime_error("error: trying to follow edges in the root");
    }

    if (get_handle_type(here) == CHAIN_HANDLE || get_handle_type(here) == SENTINEL_HANDLE) {
        assert(parent_record.get_record_handle_type() == SNARL_HANDLE ||
               parent_record.get_record_handle_type() == ROOT_HANDLE);//It could also be the root
        //If this is a chain (or a node pretending to be a chain) and it is the child of a snarl
        //Or if it is the sentinel of a snarl, then we walk through edges in the snarl
        //It can either run into another chain (or node) or the boundary node
        //TODO: What about if it is the root?


        //Get the graph handle for the end node of whatever this is, pointing in the right direction
        handle_t graph_handle;
        if (get_handle_type(here) == SENTINEL_HANDLE) {
            if ((get_connectivity(here) == START_END && !go_left) ||
                (get_connectivity(here) == START_START && go_left)) {
                graph_handle = graph->get_handle(parent_record.get_start_id(), parent_record.get_start_orientation());
            } else if ((get_connectivity(here) == END_START && !go_left) ||
                       (get_connectivity(here) == END_END && go_left)) {
                graph_handle = graph->get_handle(parent_record.get_end_id(), !parent_record.get_end_orientation());
            } else {
                //This is facing out, so don't do anything 
                return true;
            }
        } else if (get_handle_type(here) == NODE_HANDLE) {
            graph_handle = graph->get_handle(get_node_id_from_offset(get_record_offset(here)), 
                                get_connectivity(here) == START_END ? go_left : !go_left);
        } else {
            //TODO: This might not be the best way to handle orientation because it's a bit inconsistent with tips
            //Might be better to just use go_left and pretend the traversal is forward, but that might be 
            //unintuitive if you have a sentinel of a snarl that you think should be going in or out of the snarl
            //
            //If the traversal explicitly goes out the start, then we assume that it is oriented backwards
            //and go the opposite direction of go_left. Otherwise, assume that it is oriented forwards
            if (get_end_endpoint(get_connectivity(here)) == START) {
                go_left = !go_left;
            } 
            graph_handle = get_handle(get_bound(here, !go_left, false), graph);
        }
#ifdef debug_snarl_traversal
        cerr << "        traversing graph from actual node " << graph->get_id(graph_handle) << (graph->get_is_reverse(graph_handle) ? "rev" : "fd") << endl;
#endif
        graph->follow_edges(graph_handle, false, [&](const handle_t& h) {
#ifdef debug_snarl_traversal
            cerr << "  reached actual node " << graph->get_id(h) << (graph->get_is_reverse(h) ? "rev" : "fd") << endl;
#endif

            if (graph->get_id(h) == parent_record.get_start_id()) {
                //If this is the start boundary node of the parent snarl, then do this on the sentinel
                assert(graph->get_is_reverse(h) == !parent_record.get_start_orientation());
#ifdef debug_snarl_traversal
                cerr << "    -> start of parent " << endl;
#endif
                return iteratee(get_bound(get_parent(here), false, false));
            } else if (graph->get_id(h) == parent_record.get_end_id()) {
                assert(graph->get_is_reverse(h) == parent_record.get_end_orientation());
#ifdef debug_snarl_traversal
                cerr << "    -> end of parent " << endl;
#endif
                return iteratee(get_bound(get_parent(here), true, false));
            } else {
                //It is either another chain or a node, but the node needs to pretend to be a chain

                //Get the node record of the next node
                NodeRecord next_node_record(get_offset_from_node_id(graph->get_id(h)), &snarl_tree_records);

                if (next_node_record.get_parent_record_offset() == parent_record.record_offset) {
                    //If the next node's parent is also the current node's parent, then it is a node
                    //make a net_handle_t of a node pretending to be a chain
                    net_handle_t next_net = get_net_handle(next_node_record.record_offset, 
                                                           graph->get_is_reverse(h) ? END_START : START_END, 
                                                           CHAIN_HANDLE);
#ifdef debug_snarl_traversal
                cerr << "    -> actual child node " << net_handle_as_string(next_net) << endl;
#endif
                   return iteratee(next_net);
                } else {
                    //next_node_record is also the start of a chain

                    bool rev = graph->get_id(h) == next_node_record.get_is_rev_in_parent() ? false : true;
                    net_handle_t next_net = get_net_handle(next_node_record.get_parent_record_offset(), 
                                                           rev ? END_START : START_END, 
                                                           CHAIN_HANDLE);
#ifdef debug_snarl_traversal
                    assert(SnarlTreeRecord(next_node_record.get_parent_record_offset(), &snarl_tree_records).get_record_handle_type() == CHAIN_HANDLE);
                    assert(get_node_id_from_offset(next_node_record.record_offset) 
                        == SnarlTreeRecord(next_node_record.get_parent_record_offset(), &snarl_tree_records).get_start_id() || 
                        get_node_id_from_offset(next_node_record.record_offset) 
                        == SnarlTreeRecord(next_node_record.get_parent_record_offset(), &snarl_tree_records).get_end_id());
                cerr << "    -> child chain " << net_handle_as_string(next_net) << endl;
#endif
                   return iteratee(next_net);
                }
            }
            return false;
        });

        
    } else if (get_handle_type(here) == SNARL_HANDLE || get_handle_type(here) == NODE_HANDLE) {
        assert(parent_record.get_record_handle_type() == CHAIN_HANDLE);
        //If this is a snarl or node, then it is the component of a (possibly pretend) chain
        ChainRecord parent_chain(this_record.get_parent_record_offset(), &snarl_tree_records);
        net_handle_t next_net = parent_chain.get_next_child(here, go_left);
        if (next_net == here) {
            //If this is the end of the chain
            return true;
        }
#ifdef debug_snarl_traversal
        cerr << "    -> next child in the chain " << net_handle_as_string(next_net) << endl;
#endif
        return iteratee(next_net);
        
    }
    return true;
}


net_handle_t SnarlDistanceIndex::get_parent_traversal(const net_handle_t& traversal_start, const net_handle_t& traversal_end) const {
    
    net_handle_record_t start_handle_type = get_handle_type(traversal_start);
    net_handle_record_t end_handle_type = get_handle_type(traversal_end);

    if (start_handle_type == SENTINEL_HANDLE) {
        //these are the sentinels of a snarl
        //TODO: Make sure this is handling possible orientations properly
        assert(end_handle_type == SENTINEL_HANDLE);
        endpoint_t start_endpoint = get_start_endpoint(get_connectivity(traversal_start));
        endpoint_t end_endpoint = get_start_endpoint(get_connectivity(traversal_end));
        return get_net_handle(get_record_offset(get_parent(traversal_start)), 
                              endpoints_to_connectivity(start_endpoint, end_endpoint),
                              SNARL_HANDLE);
    } else {
        //These are the endpoints or tips in a chain
        SnarlTreeRecord start_record = get_snarl_tree_record(traversal_start);
        SnarlTreeRecord end_record = get_snarl_tree_record(traversal_end);
        if (start_record.get_parent_record_offset() != end_record.get_parent_record_offset()) {
            throw runtime_error("error: Looking for parent traversal of two non-siblings");
        }
        SnarlTreeRecord parent_record (start_record.get_parent_record_offset(), &snarl_tree_records);
        assert(parent_record.get_record_handle_type() == CHAIN_HANDLE);

        //Figure out what the start and end of the traversal are
        endpoint_t start_endpoint;
        if (start_handle_type == NODE_HANDLE && 
            get_node_id_from_offset(get_record_offset(traversal_start)) == parent_record.get_start_id() &&
            (get_start_endpoint(traversal_start) == START && !parent_record.get_start_orientation() ||
             get_start_endpoint(traversal_start) == END && parent_record.get_start_orientation()) ){
            //If traversal_start is a node and is also the start node oriented into the parent
            start_endpoint = START;
    
        } else if (start_handle_type == NODE_HANDLE && 
            get_node_id_from_offset(get_record_offset(traversal_start)) == parent_record.get_end_id() &&
            (get_start_endpoint(traversal_start) == START && parent_record.get_end_orientation() ||
             get_start_endpoint(traversal_start) == END && !parent_record.get_end_orientation()) ){
            //If traversal_start is a node and also the end node and oriented going into the parent
            start_endpoint = END;
    
        } else if (start_handle_type == SNARL_HANDLE) {
            start_endpoint = TIP;
        } else {
            throw runtime_error("error: trying to get an invalid traversal of a chain");
        }
    
        endpoint_t end_endpoint;
        if (end_handle_type == NODE_HANDLE && 
            get_node_id_from_offset(get_record_offset(traversal_end)) == parent_record.get_start_id() &&
            (get_start_endpoint(traversal_end) == START && parent_record.get_start_orientation() ||
             get_start_endpoint(traversal_end) == END && !parent_record.get_start_orientation())){
            //If traversal_end is a node and also the start node oriented out of the parent
            end_endpoint = START;
        } else if (end_handle_type == NODE_HANDLE && 
            get_node_id_from_offset(get_record_offset(traversal_end)) == parent_record.get_end_id() &&
            (get_start_endpoint(traversal_end) == START && !parent_record.get_end_orientation() ||
             get_start_endpoint(traversal_end) == END && parent_record.get_end_orientation()) ){
            //If traversal_end is a node and also the end node oriented out of the parent
            end_endpoint = END;
        } else if (end_handle_type == SNARL_HANDLE) {
            end_endpoint = TIP;
        } else {
            throw runtime_error("error: trying to get an invalid traversal of a chain");
        }

        if (!parent_record.has_connectivity(start_endpoint, end_endpoint)) {
            throw runtime_error("error: Trying to get parent traversal that is not connected");
        }

        return get_net_handle(parent_record.record_offset, 
                              endpoints_to_connectivity(start_endpoint, end_endpoint), 
                              CHAIN_HANDLE);
    }


}

}
