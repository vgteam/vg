//#define debug_distance_indexing
//#define debug_snarl_traversal
#define debug_distances

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {


///////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor
SnarlDistanceIndex::SnarlDistanceIndex() {}
SnarlDistanceIndex::~SnarlDistanceIndex() {}

SnarlDistanceIndex::SnarlDistanceIndex(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit) :
    snarl_size_limit(size_limit){

    //Build the temporary distance index from the graph
    TemporaryDistanceIndex temp_index(graph, snarl_finder, size_limit);

    //And fill in the permanent distance index
    vector<const TemporaryDistanceIndex*> indexes;
    indexes.emplace_back(&temp_index);
    get_snarl_tree_records(indexes, graph);
}

/*Temporary distance index for constructing the index from the graph
 */
SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryDistanceIndex(){}

SnarlDistanceIndex::TemporaryDistanceIndex::~TemporaryDistanceIndex(){}

SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryDistanceIndex(
    const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit) :
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

    //There may be components of the root that are connected to each other. Each connected component will
    //get put into a (fake) root-level snarl, but we don't know what those components will be initially, 
    //since the decomposition just puts them in the same root snarl. This is used to group the root-level
    //components into connected components that will later be used to make root snarls
    structures::UnionFind root_snarl_component_uf (0);


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
        temp_chain_record.end_node_length = graph->get_length(chain_end_handle);

        //TODO: Add root-level snarls
        if (stack.empty()) {
            //If this was the last thing on the stack, then this was a root

            //Check to see if there is anything connected to the ends of the chain
            vector<id_t> reachable_nodes;
            graph->follow_edges(graph->get_handle(temp_chain_record.start_node_id, !temp_chain_record.start_node_rev), 
                false, [&] (const handle_t& next) {
                    if (graph->get_id(next) != temp_chain_record.start_node_id &&
                        graph->get_id(next) != temp_chain_record.end_node_id) {
                        reachable_nodes.emplace_back(graph->get_id(next));
                    }
                });
            graph->follow_edges(graph->get_handle(temp_chain_record.end_node_id, temp_chain_record.end_node_rev), 
                false, [&] (const handle_t& next) {
                    if (graph->get_id(next) != temp_chain_record.start_node_id &&
                        graph->get_id(next) != temp_chain_record.end_node_id) {
                        reachable_nodes.emplace_back(graph->get_id(next));
                    }
                });
            if (reachable_nodes.size() && temp_chain_record.start_node_id != temp_chain_record.end_node_id) {
                //If we can reach anything leaving the chain (besides the chain itself), then it is part of a root snarl
                //Note that if the chain's start and end node are the same, then it will always be a single component
#ifdef debug_distance_indexing
                cerr << "                 This chain is part of the root but connects with something else in the root"<<endl;
#endif

                //Add this to the union find
                root_snarl_component_uf.resize(root_snarl_component_uf.size() + 1);
                //And remember that it's in a connected component of the root
                temp_chain_record.root_snarl_index = root_snarl_components.size();
                root_snarl_components.emplace_back(chain_index);
                for (id_t next_id : reachable_nodes) {
                    //For each node that this is connected to, check if we've already seen it and if we have, then
                    //union this chain and that node's chain
                    TemporaryNodeRecord& node_record = temp_node_records[next_id-min_node_id];
                    if (node_record.node_id != 0) {
                        //If we've already seen this node, union it with the new one
                        //If we can see it by walking out from this top-level chain, then it must also be a
                        //top-level chain (or node pretending to be a chain)
                        assert(node_record.parent.first == TEMP_CHAIN);
                        size_t other_i = temp_chain_records[node_record.parent.second].root_snarl_index;
                        assert(other_i != std::numeric_limits<size_t>::max()); 
                        root_snarl_component_uf.union_groups(other_i, temp_chain_record.root_snarl_index);
#ifdef debug_distance_indexing
                        cerr << "        Union this chain with " << temp_chain_records[node_record.parent.second].start_node_id << " " << temp_chain_records[node_record.parent.second].end_node_id << endl;
#endif
                    }
                }
            } else {
                //If this chain isn't connected to anything else, then it is a single component of the root
                temp_chain_record.parent = make_pair(TEMP_ROOT, 0);
                root_structure_count += 1;
                components.emplace_back(chain_index);
            }
        } else {
            //The last thing on the stack is the parent of this chain, which must be a snarl
            temp_chain_record.parent = stack.back();
            auto& parent_snarl_record = temp_snarl_records.at(temp_chain_record.parent.second);
            temp_chain_record.rank_in_parent = parent_snarl_record.children.size() + 2;
            parent_snarl_record.children.emplace_back(chain_index);
        }

        if (temp_chain_record.children.size() == 1 && temp_chain_record.start_node_id == temp_chain_record.end_node_id) {
            temp_chain_record.is_trivial = true;
            temp_chain_record.start_node_rev = false;
            temp_chain_record.end_node_rev = false;
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


    /*
     * We finished going through everything that exists according to the snarl decomposition, but
     * it's still missing tips, which will be discovered when filling in the snarl distances,
     * and root-level snarls, which we'll add now by combining the chain components in root_snarl_components
     * into snarls defined by root_snarl_component_uf
     * The root-level snarl is a fake snarl that doesn't exist according to the snarl decomposition, 
     * but is an extra layer that groups together components of the root that are connected
     */

    vector<vector<size_t>> root_snarl_component_indexes = root_snarl_component_uf.all_groups();
    for (vector<size_t>& root_snarl_indexes : root_snarl_component_indexes) {
#ifdef debug_distance_indexing
        cerr << "Create a new root snarl from components" << endl;
#endif
        //For each of the root snarls
        components.emplace_back(TEMP_SNARL, temp_snarl_records.size());
        temp_snarl_records.emplace_back();
        TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.back();
        temp_snarl_record.is_root_snarl = true;
        temp_snarl_record.parent = make_pair(TEMP_ROOT, 0); 

        for (size_t chain_i : root_snarl_indexes) {
            //For each chain component of this root-level snarl
            assert(root_snarl_components[chain_i].first == TEMP_CHAIN);
            TemporaryChainRecord temp_chain_record = temp_chain_records[root_snarl_components[chain_i].second];
            temp_chain_record.parent = make_pair(TEMP_SNARL, temp_snarl_records.size() - 1);
            temp_chain_record.rank_in_parent = temp_snarl_record.children.size();
            temp_chain_record.reversed_in_parent = false;

            temp_snarl_record.children.emplace_back(root_snarl_components[chain_i]);
        }
    }

    /*Now go through the decomposition again to fill in the distances
     * This traverses all chains in reverse order that we found them in, so bottom up
     * Each chain and snarl already knows its parents and children, except for single nodes
     * that are children of snarls. These nodes were not in chains will have their node 
     * records created here
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
        temp_chain_record.chain_components.emplace_back(0);


        /*First, go through each of the snarls in the chain in the forward direction and
         * fill in the distances in the snarl. Also fill in the prefix sum and backwards
         * loop vectors here
         */
        size_t curr_component = 0; //which component of the chain are we in
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

                //Fill in this snarl's distances
                populate_snarl_index(temp_snarl_record, chain_child_index, size_limit, graph);

                int64_t length;
                bool new_component = false;
                //TODO: Double check these orientations
                if (temp_snarl_record.distances.count(make_pair(make_pair(0, false), make_pair(1, false)))){
                    length = temp_snarl_record.distances.at(make_pair(make_pair(0, false), make_pair(1, false)));
                } else if (temp_snarl_record.distances.count(make_pair(make_pair(1, false), make_pair(0, false)))){
                    length = temp_snarl_record.distances.at(make_pair(make_pair(1, false), make_pair(0, false)));
                } else {
                    //The snarl is not start-end connected
                    length = std::numeric_limits<int64_t>::max();
                    //Start a new component
                    curr_component ++;
                    new_component=true;
#ifdef debug_distance_indexing
            cerr << "      This snarl is not start-end connected, starting new chain component " << endl; 
#endif
                }
                temp_snarl_record.min_length = length;



                //Get the loop distances for the snarl
                temp_snarl_record.loop_start =
                    temp_snarl_record.distances.count(make_pair(make_pair(0, false), make_pair(0, false))) 
                  ? temp_snarl_record.distances.at(make_pair(make_pair(0, false), make_pair(0, false))) 
                  : std::numeric_limits<int64_t>::max();
                temp_snarl_record.loop_end =
                    temp_snarl_record.distances.count(make_pair(make_pair(1, false), make_pair(1, false))) 
                  ? temp_snarl_record.distances.at(make_pair(make_pair(1, false), make_pair(1, false))) 
                  : std::numeric_limits<int64_t>::max();

                //And get the distance values for the end node of the snarl in the chain
                if (new_component) {
                    //If this snarl wasn't start-end connected, then we start tracking the distance vectors
                    //here
                    temp_chain_record.prefix_sum.emplace_back(0);
                    temp_chain_record.backward_loops.emplace_back(temp_snarl_record.loop_end);
                } else {
                    temp_chain_record.prefix_sum.emplace_back(sum({temp_chain_record.prefix_sum.back(), 
                        temp_snarl_record.min_length, temp_snarl_record.start_node_length}));
                    temp_chain_record.backward_loops.emplace_back(std::min(temp_snarl_record.loop_end,
                        sum({temp_chain_record.backward_loops.back() 
                        , 2 * (temp_snarl_record.start_node_length , temp_snarl_record.min_length)})));
                }
                temp_chain_record.chain_components.emplace_back(curr_component);
            }
        }
        temp_chain_record.min_length = sum({temp_chain_record.prefix_sum.back() , temp_chain_record.end_node_length}); 

        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.backward_loops.size());
        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.chain_components.size());


        /*Now that we've gone through all the snarls in the chain, fill in the forward loop vector 
         * by going through the chain in the backwards direction
         */
        temp_chain_record.forward_loops.resize(temp_chain_record.prefix_sum.size(), 
                                               std::numeric_limits<int64_t>::max());

        size_t node_i = temp_chain_record.prefix_sum.size() - 2;
        // We start at the next to last node because we need to look at this record and the next one.

        for (int j = (int)temp_chain_record.children.size() - 1 ; j >= 0 ; j--) {
            auto& child = temp_chain_record.children.at(j);
            if (child.first == TEMP_SNARL){
                TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(child.second);
                if (temp_chain_record.chain_components.at(node_i) != temp_chain_record.chain_components.at(node_i+1)){
                    //If this is a new chain component, then add the loop distance from the snarl
                    temp_chain_record.forward_loops.at(node_i) = temp_snarl_record.loop_start;
                } else if (temp_snarl_record.is_trivial) {
                    //If this is a trivial chain, then we always add the previous loop value
                    temp_chain_record.forward_loops.at(node_i) = sum({temp_chain_record.forward_loops.at(node_i+1) , 
                                    2*temp_snarl_record.end_node_length});
                } else {
                    temp_chain_record.forward_loops.at(node_i) = 
                        std::min(sum({temp_chain_record.forward_loops.at(node_i+1) , 2*temp_snarl_record.end_node_length}), 
                                 temp_snarl_record.loop_start);
                }
                node_i --;
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

    while (!all_children.empty()) {
        const pair<temp_record_t, size_t> start_index = std::move(all_children.back());
        all_children.pop_back();

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
        } else {
            assert(start_rank != 0 && start_rank != 1);
        }

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
                        node_record.is_tip = true;

                        //also update the parent
                        temp_snarl_record.node_count ++;
                        temp_snarl_record.is_trivial = false;
                        temp_snarl_record.children.emplace_back(TEMP_NODE, graph->get_id(next_handle));
                        temp_snarl_record.tippy_child_ranks.insert(node_record.rank_in_parent);

                        //TODO: Is it bad to change the list as we're walking through it?
                        all_children.emplace_back(TEMP_NODE, graph->get_id(next_handle)); 
                    }

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
                    } else {
                        assert(next_rank != 0 && next_rank != 1);
                    }
                    bool next_rev = next_index.first == TEMP_NODE || temp_chain_records[next_index.second].is_trivial 
                            ? graph->get_is_reverse(next_handle) 
                            : graph->get_id(next_handle) == temp_chain_records[next_index.second].end_node_id;

                    if (size_limit != 0 &&
                        (temp_snarl_record.node_count < size_limit ||
                         (start_rank == 0 || start_rank == 1 || next_rank == 0 || next_rank == 1))) {
                        //If we are looking at all distances or we are looking at tips or boundaries

                        //Set the distance
                        pair<size_t, bool> start = start_rank == 0 || start_rank == 1 
                            ? make_pair(start_rank, false) : make_pair(start_rank, !start_rev);
                        pair<size_t, bool> next = next_rank == 0 || next_rank == 1 
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

vector<size_t> SnarlDistanceIndex::get_snarl_tree_records(const vector<const TemporaryDistanceIndex*>& temporary_indexes, const HandleGraph* graph) {

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
    RootRecordConstructor root_record(0, total_component_count, max_node_id-min_node_id+1, min_node_id, &snarl_tree_records);
#ifdef debug_distance_indexing
    cerr << "  Root record had length " << snarl_tree_records.size() << endl;
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
                    cerr << "  Adding this chain at offset " << snarl_tree_records.size() << endl;
#endif 
                    record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records.size());

                    ChainRecordConstructor chain_record_constructor;

                    if (temp_chain_record.chain_components.back() == 0 || snarl_size_limit == 0) {
                        record_t record_type = snarl_size_limit == 0 ? CHAIN : DISTANCED_CHAIN;
                        chain_record_constructor = ChainRecordConstructor(snarl_tree_records.size(), record_type, 
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                        chain_record_constructor.set_start_end_connected();
                    } else {
                        chain_record_constructor = ChainRecordConstructor(snarl_tree_records.size(), MULTICOMPONENT_CHAIN, 
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                    }
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
                                chain_record_constructor.set_externally_start_end_connected();

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

                                //TODO: This is not really used
                                //The "rank" of the node actually points to the node in the chain, so it is the
                                //current size of the records (before adding the node to the chain)
                                node_record_constructor.set_rank_in_parent(snarl_tree_records.size());

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
                                record_to_offset.emplace(make_pair(temp_index_i, child_record_index), snarl_tree_records.size()+1);

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

                    record_t record_type = snarl_size_limit == 0 ? NODE : DISTANCED_NODE;
                    NodeRecordConstructor node_record(
                        get_offset_from_node_id(temp_node_record.node_id), record_type, &snarl_tree_records);
                    node_record.set_node_length(temp_node_record.node_length);
                    node_record.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    node_record.set_is_reversed_in_parent(false);//TODO: Make sure that this is true
                    node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);

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
                record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records.size());

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
            handle_t start_out = graph->get_handle(record.get_start_id(), !record.get_start_orientation());
            handle_t end_out = graph->get_handle(record.get_end_id(), record.get_end_orientation());
            handle_t start_in = graph->get_handle(record.get_start_id(), record.get_start_orientation());
            handle_t end_in = graph->get_handle(record.get_end_id(), !record.get_end_orientation());

            SnarlTreeRecordConstructor record_constructor(record_to_offset[make_pair(temp_index_i, component_index)], 
                                                          &snarl_tree_records);

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
            

#ifdef debug_distance_indexing
            cerr << temp_index->structure_start_end_as_string(component_index) << endl;
            assert(record.get_parent_record_offset() == 0);
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
    } else if (parent_type == ROOT_HANDLE) {
        //The parent could be a root snarl, in which case we want to actually return the root, not the 
        //snarl pretending to be the root
        return get_net_handle(0, parent_connectivity);
    }

    return get_net_handle(parent_pointer, parent_connectivity);
}

net_handle_t SnarlDistanceIndex::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    if (get_handle_type(snarl) == CHAIN_HANDLE) {
        id_t id = get_end ? SnarlTreeRecord(snarl, &snarl_tree_records).get_end_id() 
                          : SnarlTreeRecord(snarl, &snarl_tree_records).get_start_id();
        bool rev_in_parent = NodeRecord(get_offset_from_node_id(id), &snarl_tree_records).get_is_reversed_in_parent();
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
    return get_net_handle(get_record_offset(net), connectivity, get_handle_type(net));
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
        //TODO: Double check that this is the right way to handle this
        //If this is a root-level chain or node
        if ((ends_at(here) == END && !go_left) || (ends_at(here) == START && go_left)) {
            //Follow edges leaving the root structure at the end
            if (this_record.is_externally_start_end_connected()) {
                //Follow edge from end to start
                if (!iteratee(get_net_handle(get_record_offset(here), START_END, get_handle_type(here)))) {
                    return false;
                }
            }
            if (this_record.is_externally_end_end_connected()) {
                //Follow edge from end back to the end
                if (!iteratee(get_net_handle(get_record_offset(here), END_START, get_handle_type(here)))) {
                    return false;
                }
            }
        } else {
            //Follow edges leaving the root structure at the end
            if (this_record.is_externally_start_end_connected()) {
                //Follow edge from start to end
                if (!iteratee(get_net_handle(get_record_offset(here), END_START, get_handle_type(here)))) {
                    return false;
                }
            }
            if (this_record.is_externally_end_end_connected()) {
                //Follow edge from the start back to the start
                if (!iteratee(get_net_handle(get_record_offset(here), START_END, get_handle_type(here)))) {
                    return false;
                }
            }

        }
        return true;

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

                    bool rev = graph->get_id(h) == next_node_record.get_is_reversed_in_parent() ? false : true;
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

int64_t SnarlDistanceIndex::distance_in_parent(const net_handle_t& parent, 
        const net_handle_t& child1, const net_handle_t& child2) const {
#ifdef debug_distances
    cerr << "    Finding the distance in parent " << net_handle_as_string(parent) << " between " << endl 
         << "      " << net_handle_as_string(child1) << " and " <<  net_handle_as_string(child2) << endl;
#endif

    net_handle_record_t parent_handle_type = get_handle_type(parent);
    assert(canonical(parent) == canonical(get_parent(child1)));
    assert(canonical(parent) == canonical(get_parent(child2)));

    if (parent_handle_type == ROOT_HANDLE) {
        //If the parent is the root, then the children must be in the same root snarl for them to be
        //connected
        RootRecord root_record(parent, &snarl_tree_records);
        size_t parent_record_offset1 = SnarlTreeRecord(child1, &snarl_tree_records).get_parent_record_offset();
        size_t parent_record_offset2 = SnarlTreeRecord(child2, &snarl_tree_records).get_parent_record_offset();
        if (parent_record_offset1 != parent_record_offset2) {
            //If the children are in different connected components
//#ifdef debug_distances
//            cerr << "            => " << std::numeric_limits<int64_t>::max() << endl;
//#endif
            return std::numeric_limits<int64_t>::max();
        } else if (SnarlTreeRecord(parent_record_offset1, &snarl_tree_records).get_record_type() 
                    != DISTANCED_ROOT_SNARL){
            //If they are in the same connected component, but it is not a root snarl
//#ifdef debug_distances
//            cerr << "            => " << std::numeric_limits<int64_t>::max() << endl;
//#endif
            return std::numeric_limits<int64_t>::max();
        } else {
            //They are in the same root snarl, so find the distance between them
            //
            SnarlRecord snarl_record(parent_record_offset1, &snarl_tree_records);
            SnarlTreeRecord child_record1 (child1, &snarl_tree_records);
            SnarlTreeRecord child_record2 (child2, &snarl_tree_records);
            size_t rank1 = child_record1.get_rank_in_parent();
            size_t rank2 = child_record2.get_rank_in_parent();

//#ifdef debug_distances
//            cerr << "            => " << snarl_record.get_distance(rank1, ends_at(child1) == END, rank2, ends_at(child2) == END);
//#endif
            //TODO: Double check orientations
            return snarl_record.get_distance(rank1, ends_at(child1) == END, rank2, ends_at(child2) == END);
        }


    } else if (parent_handle_type == CHAIN_HANDLE) {
        ChainRecord chain_record(get_parent(child1), &snarl_tree_records);
        assert(is_node(child1) || is_snarl(child1));
        assert(is_node(child2) || is_snarl(child2));

        if (is_node(child1) && is_node(child2)) {
            //The distance between two nodes
            NodeRecord child_record1 (child1, &snarl_tree_records);
            NodeRecord child_record2 (child2, &snarl_tree_records);

            bool go_left1 = child_record1.get_is_reversed_in_parent();
            if (ends_at(child1) == START){
                go_left1 = !go_left1;
            }
            bool go_left2 = child_record2.get_is_reversed_in_parent();
            if (ends_at(child2) == START) {
                go_left2 = !go_left2;
            }
//#ifdef debug_distances
//            cerr << "            => " << chain_record.get_distance(
//                make_tuple(child_record1.get_rank_in_parent(), go_left1, child_record1.get_node_length()), 
//                make_tuple(child_record2.get_rank_in_parent(), go_left2, child_record2.get_node_length())) << endl;
//#endif

            return chain_record.get_distance(
                make_tuple(child_record1.get_rank_in_parent(), go_left1, child_record1.get_node_length()), 
                make_tuple(child_record2.get_rank_in_parent(), go_left2, child_record2.get_node_length()));
        } else if (is_node(child1) != is_node(child2)) {
            //If one of them is a node and one is a snarl
            //It doesn't matter which is which, since we're looking at the distance pointing towards each other
            NodeRecord node_record (is_node(child1) ? child1 : child2, &snarl_tree_records);
            bool go_left_node = node_record.get_is_reversed_in_parent();
            if (ends_at(is_node(child1) ? child1 : child2) == START){
                go_left_node = !go_left_node;
            }

            SnarlRecord snarl_record (is_snarl(child1) ? child1 : child2, &snarl_tree_records);
            id_t snarl_node = ends_at(is_snarl(child1) ? child1 : child2) == START 
                             ? snarl_record.get_start_id() : snarl_record.get_end_id();
            bool go_left_snarl = ends_at(is_snarl(child1) ? child1 : child2) == START ;
                                
            if (node_record.get_node_id() == snarl_node && go_left_node != go_left_snarl) {
                //If the node is the boundary of the snarl facing in
//#ifdef debug_distances
//                cerr << "        => 0" << endl;
//#endif
                return 0;
            } else {
                //Otherwise, find the actual distance from the node to the correct boundary of the snarl,
                //and add the length of the boundary of the snarl, since it is not included in the chain distance
                NodeRecord snarl_node_record (get_offset_from_node_id(snarl_node), &snarl_tree_records);
//#ifdef debug_distances
//            cerr << "            => " << sum({chain_record.get_distance(
//                               make_tuple(node_record.get_rank_in_parent(), go_left_node, node_record.get_node_length()), 
//                               make_tuple(snarl_node_record.get_rank_in_parent(), go_left_snarl, snarl_node_record.get_node_length())) 
//                           , snarl_node_record.get_node_length()}) << endl;
//#endif
                return sum({chain_record.get_distance(
                               make_tuple(node_record.get_rank_in_parent(), go_left_node, node_record.get_node_length()), 
                               make_tuple(snarl_node_record.get_rank_in_parent(), go_left_snarl, snarl_node_record.get_node_length())) 
                           , snarl_node_record.get_node_length()});
            }
        } else {
            assert(is_snarl(child1));
            assert(is_snarl(child2));
            SnarlRecord child_record1 (child1, &snarl_tree_records);
            SnarlRecord child_record2 (child2, &snarl_tree_records);

            id_t node_id1 = ends_at(child1) == START ? child_record1.get_start_id() : child_record1.get_end_id();
            id_t node_id2 = ends_at(child2) == START ? child_record2.get_start_id() : child_record2.get_end_id();
            NodeRecord node_record1 (get_offset_from_node_id(node_id1), &snarl_tree_records);
            NodeRecord node_record2 (get_offset_from_node_id(node_id2), &snarl_tree_records); 

            bool go_left1 = ends_at(child1) == START;
            bool go_left2 = ends_at(child2) == START;

            if (node_id1 == node_id2 && go_left1 != go_left2) {
                //If the snarls are adjacent (and not the same snarl)
#ifdef debug_distances
            cerr << "            => " << node_record1.get_node_length() << endl;
#endif
                return node_record1.get_node_length();
            } else {
#ifdef debug_distances
            cerr << "            => " << sum({chain_record.get_distance(
                            make_tuple(node_record1.get_rank_in_parent(), go_left1, node_record1.get_node_length()), 
                            make_tuple(node_record2.get_rank_in_parent(), go_left2, node_record2.get_node_length())) 
                    , node_record1.get_node_length() , node_record2.get_node_length()}) << endl;
#endif
                return sum({chain_record.get_distance(
                            make_tuple(node_record1.get_rank_in_parent(), go_left1, node_record1.get_node_length()), 
                            make_tuple(node_record2.get_rank_in_parent(), go_left2, node_record2.get_node_length())) 
                    , node_record1.get_node_length() , node_record2.get_node_length()});
            }
        }

    } else if (parent_handle_type == SNARL_HANDLE) {
        SnarlRecord snarl_record(get_parent(child1), &snarl_tree_records);
        SnarlTreeRecord child_record1 (child1, &snarl_tree_records);
        SnarlTreeRecord child_record2 (child2, &snarl_tree_records);
        size_t rank1, rank2; bool rev1, rev2;
        if (is_sentinel(child1)) {
            rank1 = starts_at(child1) == START ? 0 : 1;
            rev1 = false;
        } else {
            rank1 = child_record1.get_rank_in_parent();
            rev1 = ends_at(child1) == END;
        }
        if (is_sentinel(child2)) {
            rank2 = starts_at(child2) == START ? 0 : 1;
            rev2 = false;
        } else {
            rank2 = child_record2.get_rank_in_parent();
            rev2 = ends_at(child2) == END;
        }
        if ((is_sentinel(child1) && starts_at(child1) == ends_at(child1)) ||
            (is_sentinel(child2) && starts_at(child2) == ends_at(child2)) ) {
            //If this is a sentinel pointing out of the snarl
//#ifdef debug_distances
//            cerr << "            => " << std::numeric_limits<int64_t>::max() << endl;
//#endif
            return std::numeric_limits<int64_t>::max();
        }

//#ifdef debug_distances
//            cerr << "             between ranks " << rank1 << " " << rev1 << " " << rank2 << " " << rev2 << endl;
//            cerr << "            => " << snarl_record.get_distance(rank1, rev1, rank2, rev2) << endl;
//#endif
        return snarl_record.get_distance(rank1, rev1, rank2, rev2);
    } else {
        throw runtime_error("error: Trying to find distance in the wrong type of handle");
    }
}

pair<net_handle_t, bool> SnarlDistanceIndex::lowest_common_ancestor(const net_handle_t& net1, const net_handle_t& net2) const {
    net_handle_t parent1 = net1;
    net_handle_t parent2 = net2;

    std::unordered_set<net_handle_t> net1_ancestors;
    while (!is_root(parent1)){
        net1_ancestors.insert(canonical(parent1));
        parent1 = canonical(get_parent(parent1));
        cerr << "    " << net_handle_as_string(parent1) << endl;
    }

    while (net1_ancestors.count(canonical(parent2)) == 0 && !is_root(parent2)){
        //Go up until the parent2 matches something in the ancestors of net1
        //This loop will end because everything is in the same root eventually
        parent2 = canonical(get_parent(parent2));
        cerr << "    " << net_handle_as_string(parent2) << endl;
    }

    bool is_connected = true;
    if (is_root(parent2) && is_root(parent1)){
        size_t parent_record_offset1 = SnarlTreeRecord(parent1, &snarl_tree_records).get_parent_record_offset();
        size_t parent_record_offset2 = SnarlTreeRecord(parent2, &snarl_tree_records).get_parent_record_offset();
        if (parent_record_offset1 != parent_record_offset2) {
            is_connected = false;
        }
    }
    return make_pair(parent2, is_connected);
}

int64_t SnarlDistanceIndex::minimum_distance(pos_t pos1, pos_t pos2, bool unoriented_distance) const {


#ifdef debug_distances
        cerr << endl;
        cerr << "Find the minimum distance between " << pos1 << " and " << pos2 << endl;
#endif


    /*Helper function to walk up the snarl tree
     * Given a net handle, its parent,  and the distances to the start and end of the handle, 
     * update the distances to reach the ends of the parent and update the handle and its parent
     * If the parent is a chain, then the new distances include the boundary nodes of the chain.
     * If it is a snarl, it does not*/
    //TODO: This should really be an actual function and it doesn't consider the lengths of the 
    //boundary nodes. I think chains might need to know the lengths of their boundary nodes,
    //or it could find it from the snarl. Really, since the snarl is storing it anyway, I should
    //probaby rearrange things so that either the snarl or the chain can access boundary node lengths
    auto update_distances = [&](net_handle_t& net, net_handle_t& parent, int64_t& dist_start, int64_t& dist_end) {
#ifdef debug_distances
        cerr << "     at node " << net_handle_as_string(net) << " at parent " << net_handle_as_string(parent) << endl;
#endif

        if (is_trivial_chain(parent)) {
            //Don't update distances for the trivial chain
            return;
        }


        net_handle_t start_bound = get_bound(parent, false, true);
        net_handle_t end_bound = get_bound(parent, true, true);

        //The lengths of the start and end nodes of net
        //This is only needed if net is a snarl, since the boundary nodes are not technically part of the snarl
        int64_t start_node_length = 0;
        int64_t end_node_length = 0;

        //Get the distances from the bounds of the parent to the node we're looking at
        int64_t distance_start_start = distance_in_parent(parent, start_bound, flip(net));
        int64_t distance_start_end = distance_in_parent(parent, start_bound, net);
        int64_t distance_end_start = distance_in_parent(parent, end_bound, flip(net));
        int64_t distance_end_end = distance_in_parent(parent, end_bound, net);

        int64_t distance_start = dist_start;
        int64_t distance_end = dist_end; 

        int64_t start_length = is_chain(parent) ? node_length(start_bound) : 0;
        int64_t end_length = is_chain(parent) ? node_length(end_bound) : 0;

        dist_start = sum({std::min( sum({distance_start_start, distance_start}), 
                                    sum({distance_start_end , distance_end})) , 
                           start_length});
        dist_end = sum({std::min(sum({distance_end_start , distance_start}), 
                            sum({distance_end_end , distance_end})) ,
                       end_length});
#ifdef debug_distances
        cerr << "        ...new distances to start and end: " << dist_start << " " << dist_end << endl;
#endif
        return;
    };

    /*
     * Get net handles for the two nodes and the distances from each position to the ends of the handles
     * TODO: net2 is pointing in the opposite direction of the position. The final 
     * distance will be between the two nets pointing towards each other
     */
    net_handle_t net1 = get_net_handle(get_offset_from_node_id(get_id(pos1)), START_END);
    net_handle_t net2 = get_net_handle(get_offset_from_node_id(get_id(pos2)), START_END);
    pair<net_handle_t, bool> lowest_ancestor = lowest_common_ancestor(net1, net2);
    if (!lowest_ancestor.second) {
        //If these are not in the same connected component
#ifdef debug_distances
        cerr << "These are in different connected components" << endl;
#endif
        return std::numeric_limits<int64_t>::max();
    }

    //The lowest common ancestor of the two positions
    net_handle_t common_ancestor = std::move(lowest_ancestor.first);

#ifdef debug_distances
        cerr << "Found the lowest common ancestor " << net_handle_as_string(common_ancestor) << endl;
#endif
    //These are the distances to the ends of the node, including the position
    int64_t distance_to_start1 = get_is_rev(pos1) ? node_length(net1) - get_offset(pos1) : get_offset(pos1) + 1;
    int64_t distance_to_end1 = get_is_rev(pos1) ? get_offset(pos1) + 1 : node_length(net1) - get_offset(pos1);
    int64_t distance_to_start2 = get_is_rev(pos2) ? node_length(net2) - get_offset(pos2) : get_offset(pos2) + 1;
    int64_t distance_to_end2 = get_is_rev(pos2) ? get_offset(pos2) + 1 : node_length(net2) - get_offset(pos2);

    if (!unoriented_distance) {
        //If we care about the oriented distance, one of the distances will be infinite
        if (get_is_rev(pos1)) {
            distance_to_end1 = std::numeric_limits<int64_t>::max();
        } else {
            distance_to_start1 = std::numeric_limits<int64_t>::max();
        }
        if (get_is_rev(pos2)) {
            distance_to_start2 = std::numeric_limits<int64_t>::max();
        } else {
            distance_to_end2 = std::numeric_limits<int64_t>::max();
        }
    }

#ifdef debug_distances
        cerr << "Starting with distances " << distance_to_start1 << " " << distance_to_end1 << " and " << distance_to_start2 << " " << distance_to_end2 << endl;
#endif

    int64_t minimum_distance = std::numeric_limits<int64_t>::max();


    /*
     * Walk up the snarl tree until net1 and net2 are children of the lowest common ancestor
     * Keep track of the distances to the ends of the net handles as we go
     */
 
    if (canonical(net1) == canonical(net2)){
        if (get_connectivity(net1) != get_connectivity(net1) && 
            sum({distance_to_start1 , distance_to_start2}) > node_length(net1)) {
            //If the positions are on the same node and are pointing towards each other, then
            //check the distance between them in the node
            minimum_distance = minus(sum({distance_to_start1 , distance_to_start2}), node_length(net1));
        }
        common_ancestor = get_parent(net1);
    } else {

        //Get the distance from position 1 up to the ends of a child of the common ancestor
#ifdef debug_distances
        cerr << "Reaching the children of the lowest common ancestor for pos1..." << endl;
#endif   
        while (canonical(get_parent(net1)) != canonical(common_ancestor)) {
            net_handle_t parent = get_parent(net1);
            update_distances(net1, parent, distance_to_start1, distance_to_end1);
            net1 = parent;
        }
#ifdef debug_distances
        cerr << "Reached node " << net_handle_as_string(net1) << " for pos1" << endl;
        cerr << "   with distances to ends " << distance_to_start1 << " and " << distance_to_end1 << endl;
        cerr << "Reaching the children of the lowest common ancestor for pos2..." << endl;
#endif   
        //And the same for position 2
        while (canonical(get_parent(net2)) != canonical(common_ancestor)) {
            net_handle_t parent = get_parent(net2);
            update_distances(net2, parent, distance_to_start2, distance_to_end2);
            net2 = parent;
        }
#ifdef debug_distances
        cerr << "Reached node " << net_handle_as_string(net2) << " for pos2" << endl;
        cerr << "   with distances to ends " << distance_to_start2 << " and " << distance_to_end2 << endl;
#endif
    }
    net1 = canonical(net1);
    net2 = canonical(net2);


    /* 
     * common_ancestor is now the lowest common ancestor of both net handles, and 
     * net1 and net2 are both children of common_ancestor
     * Walk up to the root and check for distances between the positions within each
     * ancestor
     */
    while (!is_root(net1)){
        //TODO: Actually checking distance in a chain is between the nodes, not the snarl
        //and neither include the lengths of the nodes

        //Find the minimum distance between the two children (net1 and net2)
        int64_t distance_start_start = distance_in_parent(common_ancestor, flip(net1), flip(net2));
        int64_t distance_start_end = distance_in_parent(common_ancestor, flip(net1), net2);
        int64_t distance_end_start = distance_in_parent(common_ancestor, net1, flip(net2));
        int64_t distance_end_end = distance_in_parent(common_ancestor, net1, net2);
        cerr << "Distances: " << distance_start_start << " " << distance_start_end << " " << distance_end_start << " " << distance_end_end << endl;

        //And add those to the distances we've found to get the minimum distance between the positions
        minimum_distance = std::min(minimum_distance, 
                           std::min(sum({distance_start_start , distance_to_start1 , distance_to_start2}),
                           std::min(sum({distance_start_end , distance_to_start1 , distance_to_end2}),
                           std::min(sum({distance_end_start , distance_to_end1 , distance_to_start2}),
                                    sum({distance_end_end , distance_to_end1 , distance_to_end2})))));

#ifdef debug_distances
        cerr << "At common ancestor " << net_handle_as_string(common_ancestor) <<  endl;
        cerr << "  best distance is " << minimum_distance << endl;
#endif
        if (!is_root(common_ancestor)) {
            //Update the distances to reach the ends of the common ancestor
            update_distances(net1, common_ancestor, distance_to_start1, distance_to_end1);
            update_distances(net2, common_ancestor, distance_to_start2, distance_to_end2);

            //Update which net handles we're looking at
            net1 = common_ancestor;
            net2 = common_ancestor;
            common_ancestor = get_parent(common_ancestor);
        } else {
            //Just update this one to break out of the loop
            net1 = common_ancestor;
        }

    }

    //minimum distance currently includes both positions
    return minimum_distance == std::numeric_limits<int64_t>::max() ? std::numeric_limits<int64_t>::max() : minimum_distance-1;



}

int64_t SnarlDistanceIndex::node_length(const net_handle_t& net) const {
    assert(is_node(net));
    if (is_node(net)) {
        return NodeRecord(net, &snarl_tree_records).get_node_length();
    } else if (is_sentinel(net)) {
        //If this is the sentinel of a snarl, return the node length
        //TODO: It might be better to store the node lengths in the chains
        //This would use more memory, but it would mean a faster lookup since you wouldn't have
        //to go back down to the level of nodes. But at some point you would probably have 
        //encountered those nodes anyway, so it might not be that expensive to do it this way
        SnarlRecord snarl_record(net, &snarl_tree_records);
        NodeRecord node_record;
        if (get_start_endpoint(net) == START) {
            node_record = NodeRecord(get_offset_from_node_id(snarl_record.get_start_id()), &snarl_tree_records);
        } else {
            node_record = NodeRecord(get_offset_from_node_id(snarl_record.get_end_id()), &snarl_tree_records);
        }
        return node_record.get_node_length();
    } else {
        throw runtime_error("error: Looking for the node length of a non-node net_handle_t");
    }

}

}
