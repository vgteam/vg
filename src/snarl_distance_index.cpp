//#define debug_distance_indexing
//#define debug_snarl_traversal
//#define debug_distances

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {


///////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor
SnarlDistanceIndex::SnarlDistanceIndex() {
    snarl_tree_records.width(64);
}
SnarlDistanceIndex::~SnarlDistanceIndex() {
    snarl_tree_records.width(64);
}

SnarlDistanceIndex::SnarlDistanceIndex(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit) :
    snarl_size_limit(size_limit){
    snarl_tree_records.width(64);

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

        //Check if this is a trivial chain that represents a node
        if (temp_chain_record.children.size() == 1 && temp_chain_record.start_node_id == temp_chain_record.end_node_id) {
            temp_chain_record.is_trivial = true;
            temp_chain_record.start_node_rev = false;
            temp_chain_record.end_node_rev = false;
        }

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
            if (reachable_nodes.size() && (temp_chain_record.is_trivial || temp_chain_record.start_node_id != temp_chain_record.end_node_id)) {
                //If we can reach anything leaving the chain (besides the chain itself), then it is part of a root snarl
                //Note that if the chain's start and end node are the same, then it will always be a single component
#ifdef debug_distance_indexing
                cerr << "                 This chain is part of the root but connects with something else in the root"<<endl;
#endif
                bool new_component = true;

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
                    } else {
                        new_component = false;
                    }
                }
                if (new_component) {
                    root_structure_count += 1;
                }
            } else {
                //If this chain isn't connected to anything else, then it is a single component of the root
                temp_chain_record.parent = make_pair(TEMP_ROOT, 0);
                components.emplace_back(chain_index);
                root_structure_count += 1;
            }
        } else {
            //The last thing on the stack is the parent of this chain, which must be a snarl
            temp_chain_record.parent = stack.back();
            auto& parent_snarl_record = temp_snarl_records.at(temp_chain_record.parent.second);
            temp_chain_record.rank_in_parent = parent_snarl_record.children.size() + 2;
            parent_snarl_record.children.emplace_back(chain_index);
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
            TemporaryChainRecord& temp_chain_record = temp_chain_records[root_snarl_components[chain_i].second];
            temp_chain_record.parent = make_pair(TEMP_SNARL, temp_snarl_records.size() - 1);
            temp_chain_record.rank_in_parent = temp_snarl_record.children.size();
            temp_chain_record.reversed_in_parent = false;

            temp_snarl_record.children.emplace_back(root_snarl_components[chain_i]);
        }
        temp_snarl_record.node_count = temp_snarl_record.children.size();
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
        for (size_t chain_child_i = 0 ; chain_child_i < temp_chain_record.children.size() ; chain_child_i++ ){ 
            const pair<temp_record_t, size_t>& chain_child_index = temp_chain_record.children[chain_child_i];
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

                    //If this is the second component of the multicomponent chain, then remember the minimum length
                    if (curr_component == 1) {
                        temp_chain_record.min_length = temp_chain_record.prefix_sum.back();
                    }
                    temp_chain_record.prefix_sum.emplace_back(0);
                    temp_chain_record.backward_loops.emplace_back(temp_snarl_record.loop_end);
                } else {
                    temp_chain_record.prefix_sum.emplace_back(sum({temp_chain_record.prefix_sum.back(), 
                        temp_snarl_record.min_length, temp_snarl_record.start_node_length}));
                    temp_chain_record.backward_loops.emplace_back(std::min(temp_snarl_record.loop_end,
                        sum({temp_chain_record.backward_loops.back() 
                        , 2 * (temp_snarl_record.start_node_length + temp_snarl_record.min_length)})));
                }
                temp_chain_record.chain_components.emplace_back(curr_component);
                if (chain_child_i == temp_chain_record.children.size() - 2 && temp_snarl_record.min_length == std::numeric_limits<int64_t>::max()) {
                    temp_chain_record.loopable = false;
                }
            }
        } //Finished walking through chain
        if (temp_chain_record.start_node_id == temp_chain_record.end_node_id && temp_chain_record.chain_components.back() != 0) {
            //If this is a looping, multicomponent chain, the start/end node could end up in separate chain components
            //despite being the same node.
            //Since the first component will always be 0, set the first node's component to be whatever the last 
            //component was
            temp_chain_record.chain_components[0] = temp_chain_record.chain_components.back();

        }
        temp_chain_record.min_length = !temp_chain_record.is_trivial && temp_chain_record.start_node_id == temp_chain_record.end_node_id 
                        ? sum({temp_chain_record.prefix_sum.back(), temp_chain_record.min_length})
                        : sum({temp_chain_record.prefix_sum.back() , temp_chain_record.end_node_length});

        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.backward_loops.size());
        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.chain_components.size());


        /*Now that we've gone through all the snarls in the chain, fill in the forward loop vector 
         * by going through the chain in the backwards direction
         */
        temp_chain_record.forward_loops.resize(temp_chain_record.prefix_sum.size(), 
                                               std::numeric_limits<int64_t>::max());
        if (temp_chain_record.start_node_id == temp_chain_record.end_node_id && temp_chain_record.children.size() > 1) {
            //If this is a looping chain, then check the first snarl for a loop
            TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(temp_chain_record.children.at(1).second);
            temp_chain_record.forward_loops[temp_chain_record.forward_loops.size()-1] = temp_snarl_record.loop_start;


            //Also check if the reverse loop values would be improved if we went around again
            if (temp_chain_record.backward_loops.back() < temp_chain_record.backward_loops.front()) {
                temp_chain_record.backward_loops[0] = temp_chain_record.backward_loops.back();
            }
        }

        size_t node_i = temp_chain_record.prefix_sum.size() - 2;
        // We start at the next to last node because we need to look at this record and the next one.

        for (int j = (int)temp_chain_record.children.size() - 1 ; j >= 0 ; j--) {
            auto& child = temp_chain_record.children.at(j);
            if (child.first == TEMP_SNARL){
                TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(child.second);
                if (temp_chain_record.chain_components.at(node_i) != temp_chain_record.chain_components.at(node_i+1)){
                    //If this is a new chain component, then add the loop distance from the snarl
                    temp_chain_record.forward_loops.at(node_i) = temp_snarl_record.loop_start;
                } else {
                    temp_chain_record.forward_loops.at(node_i) = 
                        std::min(sum({temp_chain_record.forward_loops.at(node_i+1) , 2* temp_snarl_record.min_length,
                                      2*temp_snarl_record.end_node_length}), temp_snarl_record.loop_start);
                }
                node_i --;
            }
        }
    }

#ifdef debug_distance_indexing
    cerr << "Filling in the distances in root snarls" << endl;
#endif
    for (pair<temp_record_t, size_t>& component_index : components) {
        if (component_index.first == TEMP_SNARL) {
            TemporarySnarlRecord& temp_snarl_record = temp_snarl_records.at(component_index.second);
            populate_snarl_index(temp_snarl_record, component_index, size_limit, graph);
            temp_snarl_record.min_length = std::numeric_limits<int64_t>::max();//TODO: This is true but might be better to store it as something else so we can bit compress later
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

bdsg::MappedIntVector SnarlDistanceIndex::get_snarl_tree_records(const vector<const TemporaryDistanceIndex*>& temporary_indexes, const HandleGraph* graph) {

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
                    cerr << "            with indices " << current_record_index.first << " " << current_record_index.second << endl;
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
                cerr << "       child " << temp_index->structure_start_end_as_string(child) << endl;
                cerr << "        " << child.first << " " << child.second << endl;
                cerr << "        Add child " << net_handle_as_string(get_net_handle(record_to_offset[make_pair(temp_index_i, child)], START_END)) 
                     << "     at offset " << record_to_offset[make_pair(temp_index_i, child)] 
                     << "     to child list at offset " << snarl_tree_records.size() << endl;
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
    return get_handle_type(net) == ROOT_HANDLE;
}

bool SnarlDistanceIndex::is_snarl(const net_handle_t& net) const {
    return get_handle_type(net) == SNARL_HANDLE 
            && (SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == SNARL_HANDLE ||
                SnarlTreeRecord(net, &snarl_tree_records).get_record_type() == ROOT_SNARL ||  
                SnarlTreeRecord(net, &snarl_tree_records).get_record_type() == DISTANCED_ROOT_SNARL);
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
                                 ends_at(net) == END ? false : true);
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
        bool rev_in_parent = NodeRecord(get_offset_from_node_id(id), &snarl_tree_records).get_is_reversed_in_parent();
        if (get_end) {
            rev_in_parent = !rev_in_parent;
        }
        if (!face_in){
            rev_in_parent = !rev_in_parent;
        }
        connectivity_t connectivity = rev_in_parent ? END_START : START_END;
        if (SnarlTreeRecord(snarl, &snarl_tree_records).get_start_id() 
                == SnarlTreeRecord(snarl, &snarl_tree_records).get_end_id() &&  get_end) {
            //If this is a looping chain and we're getting the end, then the traversal of the node will be the 
            //end endpoint repeated, instead of the actual traversal
            //TODO: This might cause problems when checking traversals but I think it's fine
            connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), get_end_endpoint(connectivity));
        }
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

    } else if (get_handle_type(here) == CHAIN_HANDLE || get_handle_type(here) == SENTINEL_HANDLE) {
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
                                ends_at(here) == END ? go_left : !go_left);
        } else {
            //TODO: This might not be the best way to handle orientation because it's a bit inconsistent with tips
            //Might be better to just use go_left and pretend the traversal is forward, but that might be 
            //unintuitive if you have a sentinel of a snarl that you think should be going in or out of the snarl
            //
            //If the traversal explicitly goes out the start, then we assume that it is oriented backwards
            //and go the opposite direction of go_left. Otherwise, assume that it is oriented forwards
            if (ends_at(here) == START) {
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
        if (ends_at(here) == START) {
            go_left = !go_left;
        }
        if (SnarlTreeRecord(get_record_offset(here), &snarl_tree_records).get_is_reversed_in_parent()) {
            go_left = !go_left;
        }
        net_handle_t next_net = parent_chain.get_next_child(here, go_left);
        if (next_net == here) {
            //If this is the end of the chain
            return true;
        }
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
        const net_handle_t& child1, const net_handle_t& child2, const HandleGraph* graph) const {

    assert(canonical(parent) == canonical(get_parent(child1)));
    assert(canonical(parent) == canonical(get_parent(child2)));

    if (is_root(parent)) {
        //If the parent is the root, then the children must be in the same root snarl for them to be
        //connected
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


    } else if (is_chain(parent)) {
        ChainRecord chain_record(get_parent(child1), &snarl_tree_records);
        assert(is_node(child1) || is_snarl(child1));
        assert(is_node(child2) || is_snarl(child2));

        if (is_trivial_chain(parent)) {
            return std::numeric_limits<int64_t>::max();
        }

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
//            cerr << "Finding distances between two nodes with ranks " << child_record1.get_rank_in_parent() << " " << go_left1 << " " << child_record1.get_node_length() << " and " << endl << child_record2.get_rank_in_parent() << " " << go_left2 << " " << child_record2.get_node_length() << endl;
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
//#ifdef debug_distances
//            cerr << "            => " << node_record1.get_node_length() << endl;
//#endif

                return node_record1.get_node_length();
            } else {
//#ifdef debug_distances
//            cerr << "            => " << sum({chain_record.get_distance(
//                            make_tuple(node_record1.get_rank_in_parent(), go_left1, node_record1.get_node_length()), 
//                            make_tuple(node_record2.get_rank_in_parent(), go_left2, node_record2.get_node_length())) 
//                    , node_record1.get_node_length() , node_record2.get_node_length()}) << endl;
//#endif
                return sum({chain_record.get_distance(
                            make_tuple(node_record1.get_rank_in_parent(), go_left1, node_record1.get_node_length()), 
                            make_tuple(node_record2.get_rank_in_parent(), go_left2, node_record2.get_node_length())) 
                    , node_record1.get_node_length() , node_record2.get_node_length()});
            }
        }

    } else if (is_snarl(parent)) {
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
//
        if (snarl_record.get_record_type() == OVERSIZED_SNARL 
            && !(rank1 == 0 || rank1 == 1 || rank2 == 0 || rank2 == 1)) {
            //If this is an oversized snarl and we're looking for internal distances, then we didn't store the
            //distance and we have to find it using dijkstra's algorithm
            if (graph == nullptr) {
                cerr << "warning: trying to find the distance in an oversized snarl without a graph. Returning inf" << endl;
                return std::numeric_limits<int64_t>::max();
            }
            handle_t handle1 = is_node(child1) ? get_handle(child1, graph) : get_handle(get_bound(child1, ends_at(child1) == END, false), graph); 
            handle_t handle2 = is_node(child2) ? get_handle(child2, graph) : get_handle(get_bound(child2, ends_at(child2) == END, false), graph);
            handle2 = graph->flip(handle2);

            int64_t distance = std::numeric_limits<int64_t>::max();
            handlegraph::algorithms::dijkstra(graph, handle1, [&](const handle_t& reached, size_t dist) {
                if (reached == handle2) {
                    //TODO: Also give up if the distance is too great
                    distance = dist;
                    return false;
                }
                return true;
            }, false);
            return distance;

            
        } else {
           return snarl_record.get_distance(rank1, rev1, rank2, rev2);
        }
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
    }

    while (net1_ancestors.count(canonical(parent2)) == 0 && !is_root(parent2)){
        //Go up until the parent2 matches something in the ancestors of net1
        //This loop will end because everything is in the same root eventually
        parent2 = canonical(get_parent(parent2));
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

int64_t SnarlDistanceIndex::minimum_distance(pos_t pos1, pos_t pos2, bool unoriented_distance, const HandleGraph* graph) const {


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
        cerr << "     Updating distance from node " << net_handle_as_string(net) << " at parent " << net_handle_as_string(parent) << endl;
#endif

        if (is_trivial_chain(parent)) {
            //Don't update distances for the trivial chain
            return;
        }


        net_handle_t start_bound = get_bound(parent, false, true);
        net_handle_t end_bound = get_bound(parent, true, true);

        //The lengths of the start and end nodes of net
        //This is only needed if net is a snarl, since the boundary nodes are not technically part of the snarl
        int64_t start_length = is_chain(parent) ? node_length(start_bound) : 0;
        int64_t end_length = is_chain(parent) ? node_length(end_bound) : 0;

        //Get the distances from the bounds of the parent to the node we're looking at
        int64_t distance_start_start = start_bound == net ? -start_length : distance_in_parent(parent, start_bound, flip(net), graph);
        int64_t distance_start_end = start_bound == flip(net) ? -start_length : distance_in_parent(parent, start_bound, net, graph);
        int64_t distance_end_start = end_bound == net ? -end_length : distance_in_parent(parent, end_bound, flip(net), graph);
        int64_t distance_end_end = end_bound == flip(net) ? -end_length : distance_in_parent(parent, end_bound, net, graph);

        int64_t distance_start = dist_start;
        int64_t distance_end = dist_end; 


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
        if (sum({distance_to_end1 , distance_to_start2}) > node_length(net1) && 
            sum({distance_to_end1 , distance_to_start2}) != std::numeric_limits<int64_t>::max()) {
            //If the positions are on the same node and are pointing towards each other, then
            //check the distance between them in the node
            minimum_distance = minus(sum({distance_to_end1 , distance_to_start2}), node_length(net1));
        }
        if (sum({distance_to_start1 , distance_to_end2}) > node_length(net1) && 
            sum({distance_to_start1 , distance_to_end2}) != std::numeric_limits<int64_t>::max()) {
            minimum_distance = std::min(minus(sum({distance_to_start1 , distance_to_end2}), node_length(net1)), minimum_distance);
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
    //TODO: I'm taking this out because it should really be start-end connected, but this
    //won't do this if that traversal isn't possible. Really it should just be setting the 
    //connectivity instead
    //net1 = canonical(net1);
    //net2 = canonical(net2);


    /* 
     * common_ancestor is now the lowest common ancestor of both net handles, and 
     * net1 and net2 are both children of common_ancestor
     * Walk up to the root and check for distances between the positions within each
     * ancestor
     */
    while (!is_root(net1)){
        //TODO: Actually checking distance in a chain is between the nodes, not the snarl
        //and neither include the lengths of the nodes
#ifdef debug_distances
            cerr << "At common ancestor " << net_handle_as_string(common_ancestor) <<  endl;
            cerr << "  with distances for child 1 (" << net_handle_as_string(net1) << "): " << distance_to_start1 << " "  << distance_to_end1 << endl;
            cerr << "                     child 2 (" << net_handle_as_string(net2) << "): " << distance_to_start2 << " " <<  distance_to_end2 << endl;
#endif

        //Find the minimum distance between the two children (net1 and net2)
        int64_t distance_start_start = distance_in_parent(common_ancestor, flip(net1), flip(net2), graph);
        int64_t distance_start_end = distance_in_parent(common_ancestor, flip(net1), net2, graph);
        int64_t distance_end_start = distance_in_parent(common_ancestor, net1, flip(net2), graph);
        int64_t distance_end_end = distance_in_parent(common_ancestor, net1, net2, graph);

        //And add those to the distances we've found to get the minimum distance between the positions
        minimum_distance = std::min(minimum_distance, 
                           std::min(sum({distance_start_start , distance_to_start1 , distance_to_start2}),
                           std::min(sum({distance_start_end , distance_to_start1 , distance_to_end2}),
                           std::min(sum({distance_end_start , distance_to_end1 , distance_to_start2}),
                                    sum({distance_end_end , distance_to_end1 , distance_to_end2})))));

#ifdef debug_distances
            cerr << "    Found distances between nodes: " << distance_start_start << " " << distance_start_end << " " << distance_end_start << " " << distance_end_end << endl;
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
            //If net1 and net2 are the same, check if they have external connectivity in the root
            if ( canonical(net1) == canonical(net2) ) {
#ifdef debug_distances
                cerr << "    Checking external connectivity" << endl;
#endif
                if (has_external_connectivity(net1, START, START)) {
                    minimum_distance = std::min(minimum_distance, sum({distance_to_start1, distance_to_start2}));
                }
                if (has_external_connectivity(net1, END, END)) {
                    minimum_distance = std::min(minimum_distance, sum({distance_to_end1, distance_to_end2}));
                }
                if (has_external_connectivity(net1, START, END)) {
                    minimum_distance = std::min(minimum_distance, 
                                       std::min(sum({distance_to_start1, distance_to_end2}),
                                                sum({distance_to_end1, distance_to_start2})));
                }

                if (has_external_connectivity(net1, START, START) && has_external_connectivity(net1, END, END)) {
                    minimum_distance = std::min(minimum_distance, 
                                       sum({std::min(sum({distance_to_start1, distance_to_end2}),
                                                sum({distance_to_end1, distance_to_start2})), 
                                           minimum_length(net1)}));
                }
            }
            //Just update this one to break out of the loop
            net1 = common_ancestor;

        }

#ifdef debug_distances
            cerr << "  new common ancestor " << net_handle_as_string(common_ancestor) <<  endl;
            cerr << "  new distances are " << distance_to_start1 << " "  << distance_to_end1 << 
                    " " << distance_to_start2 << " " <<  distance_to_end2 << endl;
#endif
    }

    //minimum distance currently includes both positions
    return minimum_distance == std::numeric_limits<int64_t>::max() ? std::numeric_limits<int64_t>::max() : minimum_distance-1;



}

int64_t SnarlDistanceIndex::node_length(const net_handle_t& net) const {
    assert(is_node(net) || is_sentinel(net));
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


//TODO: This is kind of redundant with node_length 
int64_t SnarlDistanceIndex::minimum_length(const net_handle_t& net) const {
        return SnarlTreeRecord(net, &snarl_tree_records).get_min_length();
}
int64_t SnarlDistanceIndex::node_id(const net_handle_t& net) const {
    assert(is_node(net) || is_sentinel(net));
    if (is_node(net)) {
        return NodeRecord(net, &snarl_tree_records).get_node_id();
    } else if (is_sentinel(net)) {
        SnarlRecord snarl_record(net, &snarl_tree_records);
        NodeRecord node_record;
        if (get_start_endpoint(net) == START) {
            return snarl_record.get_start_id();
        } else {
            return snarl_record.get_end_id();
        }
    } else {
        throw runtime_error("error: Looking for the node length of a non-node net_handle_t");
    }

}

bool SnarlDistanceIndex::has_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    return record.has_connectivity(start, end);
}
bool SnarlDistanceIndex::has_external_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    return record.has_external_connectivity(start, end);
}
bool SnarlDistanceIndex::SnarlTreeRecord::has_connectivity(connectivity_t connectivity) const {
    if (connectivity == START_START) {
        return is_start_start_connected();
    } else if (connectivity == START_END || connectivity == END_START) {
        return is_start_end_connected();
    } else if (connectivity == START_TIP || connectivity == TIP_START) {
        return is_start_tip_connected();
    } else if (connectivity == END_END) {
        return is_end_end_connected();
    } else if (connectivity == END_TIP || connectivity == TIP_END) {
        return is_end_tip_connected();
    } else if (connectivity == TIP_TIP) {
        return is_tip_tip_connected();
    } else {
        throw runtime_error("error: Invalid connectivity");
    }
}
size_t SnarlDistanceIndex::SnarlTreeRecord::get_min_length() const {
    record_t type = get_record_type();
    if (type == DISTANCED_NODE) {
        return records->at(record_offset + NODE_LENGTH_OFFSET);
    } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return records->at(record_offset + SNARL_MIN_LENGTH_OFFSET);
    } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return records->at(record_offset + CHAIN_MIN_LENGTH_OFFSET);
    } else if (type == NODE || type == SNARL || type == CHAIN) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (type == ROOT || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the length of the root");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

size_t SnarlDistanceIndex::SnarlTreeRecord::get_max_length() const {
    record_t type = get_record_type();
    if (type == DISTANCED_NODE) {
        return records->at(record_offset + NODE_LENGTH_OFFSET);
    } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return records->at(record_offset + SNARL_MAX_LENGTH_OFFSET);
    } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return records->at(record_offset + CHAIN_MAX_LENGTH_OFFSET);
    } else if (type == NODE || type == SNARL || type == CHAIN) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }  else if (type == ROOT || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the length of the root");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

size_t SnarlDistanceIndex::SnarlTreeRecord::get_rank_in_parent() const {
    record_t type = get_record_type();
    if (type == NODE || type == DISTANCED_NODE) {
        return records->at(record_offset + NODE_RANK_OFFSET) >> 1;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
        return records->at(record_offset + SNARL_RANK_OFFSET) >> 1;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return records->at(record_offset + CHAIN_RANK_OFFSET) >> 1;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};
bool SnarlDistanceIndex::SnarlTreeRecord::get_is_reversed_in_parent() const {
    record_t type = get_record_type();
    if (type == NODE || type == DISTANCED_NODE) {
        return records->at(record_offset + NODE_RANK_OFFSET) & 1;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
        return records->at(record_offset + SNARL_RANK_OFFSET) & 1;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return records->at(record_offset + CHAIN_RANK_OFFSET) & 1;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

id_t SnarlDistanceIndex::SnarlTreeRecord::get_start_id() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        //TODO: Also not totally sure what this should do
        throw runtime_error("error: trying to get the start node of the root");
    } else if (type == NODE || type == DISTANCED_NODE) {
        //TODO: I'm not sure if I want to allow this
        //cerr << "warning: Looking for the start of a node" << endl;
        return get_id_from_offset(record_offset);
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return (records->at(record_offset + SNARL_START_NODE_OFFSET)) >> 1;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return (records->at(record_offset + CHAIN_START_NODE_OFFSET)) >> 1;
    } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the start node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}
bool SnarlDistanceIndex::SnarlTreeRecord::get_start_orientation() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        //TODO: Also not totally sure what this should do
        throw runtime_error("error: trying to get the start node of the root");
    } else if (type == NODE || type == DISTANCED_NODE) {
        //TODO: I'm not sure if I want to allow this
        //cerr << "warning: Looking for the start of a node" << endl;
        return false;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return (records->at(record_offset + SNARL_START_NODE_OFFSET) & 1);
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return (records->at(record_offset + CHAIN_START_NODE_OFFSET)) & 1;
    } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the start node of a root snarl");
    }else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}
id_t SnarlDistanceIndex::SnarlTreeRecord::get_end_id() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        //TODO: Also not totally sure what this should do
        throw runtime_error("error: trying to get the end node of the root");
    } else if (type == NODE || type == DISTANCED_NODE) {
        //TODO: I'm not sure if I want to allow this
        //cerr << "warning: Looking for the end of a node" << endl;
        //TODO: Put this in its own function? Also double check for off by ones
        //Offset of the start of the node vector
        return get_id_from_offset(record_offset);
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return (records->at(record_offset + SNARL_END_NODE_OFFSET)) >> 1;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return (records->at(record_offset + CHAIN_END_NODE_OFFSET)) >> 1;
    } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}

id_t SnarlDistanceIndex::SnarlTreeRecord::get_end_orientation() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        //TODO: Also not totally sure what this should do
        throw runtime_error("error: trying to get the end node of the root");
    } else if (type == NODE || type == DISTANCED_NODE) {
        //TODO: I'm not sure if I want to allow this
        //cerr << "warning: Looking for the end of a node" << endl;
        //TODO: Put this in its own function? Also double check for off by ones
        //Offset of the start of the node vector
        return get_id_from_offset(record_offset);
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        return (records->at(record_offset + SNARL_END_NODE_OFFSET)) & 1;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return (records->at(record_offset + CHAIN_END_NODE_OFFSET)) & 1;
    } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: trying to find the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}

size_t SnarlDistanceIndex::SnarlTreeRecord::get_offset_from_id (const id_t id) const {
    size_t node_records_offset = records->at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE;
    size_t offset = (id-records->at(MIN_NODE_ID_OFFSET)) * NODE_RECORD_SIZE;
    return node_records_offset + offset;

}
id_t SnarlDistanceIndex::SnarlTreeRecord::get_id_from_offset(size_t offset) const {
    size_t min_node_id = records->at(MIN_NODE_ID_OFFSET);
    size_t node_records_offset = records->at(COMPONENT_COUNT_OFFSET) + ROOT_RECORD_SIZE;
    return ((offset-node_records_offset) / NODE_RECORD_SIZE) + min_node_id;
}


bool SnarlDistanceIndex::SnarlTreeRecord::has_connectivity(endpoint_t start, endpoint_t end){
    return has_connectivity(endpoints_to_connectivity(start, end));
}

bool SnarlDistanceIndex::SnarlTreeRecord::has_external_connectivity(endpoint_t start, endpoint_t end) {
    if (start == START && end == START) {
        return is_externally_start_start_connected();
    } else if (start == END && end == END) {
        return is_externally_end_end_connected();
    } else if ( (start == START && end == END ) || (start == END && end == START)) {
        return is_externally_start_end_connected();
    } else {
        return false;
    }
}

size_t SnarlDistanceIndex::SnarlTreeRecord::get_parent_record_offset() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        return 0;
    } else if (type == NODE || type == DISTANCED_NODE) {
        return (records->at(record_offset + NODE_PARENT_OFFSET));
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
        return (records->at(record_offset + SNARL_PARENT_OFFSET));
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        return (records->at(record_offset + CHAIN_PARENT_OFFSET));
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

SnarlDistanceIndex::record_t SnarlDistanceIndex::SnarlTreeRecordConstructor::get_record_type() const {
    return static_cast<record_t>(records->at(record_offset) >> 9);
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_start_start_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set start_start_connected" << endl;
#endif    

    records->at(record_offset) = records->at(record_offset) | 32;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_start_end_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set start_end_connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 16;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_start_tip_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set start_tip_connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 8;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_end_end_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set end_end_connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 4;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_end_tip_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set end_tip_connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 2;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_tip_tip_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set tpi_tip_connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 1;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_externally_start_end_connected() {
#ifdef debug_indexing
    cerr << record_offset << " set externally start_end connected" << endl;
#endif

    records->at(record_offset) = records->at(record_offset) | 64;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_externally_start_start_connected() const {
#ifdef debug_indexing
    cerr << record_offset << " set externally start_start connected" << endl;
#endif
    records->at(record_offset) = records->at(record_offset) | 128;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_externally_end_end_connected() const {
#ifdef debug_indexing
    cerr << record_offset << " set externally end_end connected" << endl;
#endif
    records->at(record_offset) = records->at(record_offset) | 256;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_record_type(record_t type) {
#ifdef debug_indexing
    cerr << record_offset << " set record type to be " << type << endl;
    assert(records->at(record_offset) == 0);
#endif
    records->at(record_offset) = ((static_cast<size_t>(type) << 9) | (records->at(record_offset) & 511));
}


void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_min_length(size_t length) {
    record_t type = get_record_type();
    size_t offset;
    if (type == DISTANCED_NODE) {
        offset = record_offset + NODE_LENGTH_OFFSET;
    } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        offset = record_offset + SNARL_MIN_LENGTH_OFFSET;
    } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_MIN_LENGTH_OFFSET;
    } else if (type == NODE || type == SNARL || type == CHAIN ) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: set the length of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set min length to be " << length << endl;
    assert(records->at(offset) == 0);
#endif

    records->at(offset) = length;
};
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_max_length(size_t length) {
    record_t type = get_record_type();
    size_t offset;
    if (type == DISTANCED_NODE) {
        offset = record_offset + NODE_LENGTH_OFFSET;
    } else if (type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        offset = record_offset + SNARL_MAX_LENGTH_OFFSET;
    } else if (type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_MAX_LENGTH_OFFSET;
    } else if (type == DISTANCED_NODE || type == SNARL || type == CHAIN) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: set the length of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set max length to be " << length << endl;
    assert(records->at(offset) == 0);
#endif

    records->at(offset) = length;
};

void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_rank_in_parent(size_t rank) {
    record_t type = get_record_type();
    size_t offset;
    if (type == NODE || type == DISTANCED_NODE) {
        offset = record_offset + NODE_RANK_OFFSET;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
        offset = record_offset + SNARL_RANK_OFFSET;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_RANK_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set rank in parent to be " << rank << endl;
    assert(records->at(offset) >> 1 == 0);
#endif

    bool rev = records->at(offset) & 1;
    records->at(offset) = (rank << 1) | rev;
};
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_is_reversed_in_parent(bool rev) {
    record_t type = get_record_type();
    size_t offset;
    if (type == NODE || type == DISTANCED_NODE) {
        offset = record_offset + NODE_RANK_OFFSET;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
        offset = record_offset + SNARL_RANK_OFFSET;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_RANK_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set rev in parent to be " << rev << endl;
#endif

    records->at(offset) =  ((records->at(offset)>>1)<<1) | rev;
};
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_parent_record_offset(size_t pointer){
    record_t type = get_record_type();
    size_t offset;
    if (type == NODE || type == DISTANCED_NODE) {
        offset = record_offset + NODE_PARENT_OFFSET;
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL
            || type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL)  {
#ifdef debug_indexing
        if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
            assert(pointer == 0);
        }
#endif

        offset = record_offset + SNARL_PARENT_OFFSET;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_PARENT_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set parent offset to be " << pointer << endl;
    assert(records->at(offset) == 0);
#endif

    records->at(offset) = pointer;
};
//Rev is true if the node is traversed backwards to enter the snarl
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_start_node(id_t id, bool rev) {
    record_t type = get_record_type();
    size_t offset;
    if (type == ROOT || type == NODE || type == DISTANCED_NODE) {
        throw runtime_error("error: trying to set the node id of a node or root");
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        offset = record_offset + SNARL_START_NODE_OFFSET;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_START_NODE_OFFSET;
    }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: set the start node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set start node to be " << id << " facing " << (rev ? "rev" : "fd") << endl;
    assert(records->at(offset) == 0);
#endif

    records->at(offset) = (id << 1) | rev;
}
void SnarlDistanceIndex::SnarlTreeRecordConstructor::set_end_node(id_t id, bool rev) const {
    record_t type = get_record_type();
    size_t offset;
    if (type == ROOT || type == NODE || type == DISTANCED_NODE) {
        throw runtime_error("error: trying to set the node id of a node or root");
    } else if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL)  {
        offset = record_offset + SNARL_END_NODE_OFFSET;
    } else if (type == CHAIN || type == DISTANCED_CHAIN || type == MULTICOMPONENT_CHAIN)  {
        offset = record_offset + CHAIN_END_NODE_OFFSET;
    }  else if (type == ROOT_SNARL || type == DISTANCED_ROOT_SNARL) {
        throw runtime_error("error: set the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
#ifdef debug_indexing
    cerr << offset << " set end node to be " << id << " facing " << (rev ? "rev" : "fd") << endl;
    assert(records->at(offset) == 0);
#endif

    records->at(offset) = (id << 1) | rev;
}



bool SnarlDistanceIndex::RootRecord::for_each_child(const std::function<bool(const handlegraph::net_handle_t&)>& iteratee) const {
    size_t connected_component_count = get_connected_component_count();
    for (size_t i = 0 ; i < connected_component_count ; i++) {

        size_t child_offset = records->at(record_offset + ROOT_RECORD_SIZE + i);
        net_handle_record_t type = SnarlTreeRecord(child_offset, records).get_record_handle_type();
        record_t record_type = SnarlTreeRecord(child_offset, records).get_record_type();


        if (record_type == ROOT_SNARL || record_type == DISTANCED_ROOT_SNARL) {
            //This is a bunch of root components that are connected, so go through each
            SnarlRecord snarl_record(child_offset, records);
            if (! snarl_record.for_each_child(iteratee)) {
                return false;
            }
        } else {
            //Otherwise, it is a separate connected component
            net_handle_t child_handle =  get_net_handle(child_offset, START_END, type);
            if (!iteratee(child_handle)) {
                return false;
            }
        }
    }
    return true;
}

void SnarlDistanceIndex::RootRecordConstructor::set_connected_component_count(size_t connected_component_count) {
#ifdef debug_indexing
    cerr << RootRecord::record_offset+COMPONENT_COUNT_OFFSET << " set connected component to be " << connected_component_count << endl;
    assert(SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+COMPONENT_COUNT_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+COMPONENT_COUNT_OFFSET)=connected_component_count;
}
void SnarlDistanceIndex::RootRecordConstructor::set_node_count(size_t node_count) {
#ifdef debug_indexing
    cerr << RootRecord::record_offset+NODE_COUNT_OFFSET << " set node count to be " << node_count << endl;
    assert(RootRecord::records->at(RootRecord::record_offset+NODE_COUNT_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+NODE_COUNT_OFFSET)=node_count;
}
void SnarlDistanceIndex::RootRecordConstructor::set_min_node_id(id_t node_id) {
#ifdef debug_indexing
    cerr << RootRecord::record_offset+MIN_NODE_ID_OFFSET << " set min node id to be " << node_id << endl;
    assert(RootRecord::records->at(RootRecord::record_offset+MIN_NODE_ID_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+MIN_NODE_ID_OFFSET)=node_id;
}
void SnarlDistanceIndex::RootRecordConstructor::add_component(size_t index, size_t offset) {
#ifdef debug_indexing
    cerr << RootRecord::record_offset+ROOT_RECORD_SIZE+index << " set new component " << offset << endl;
    assert(SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+ROOT_RECORD_SIZE+index) == 0);
    assert(index < get_connected_component_count());
#endif

    SnarlTreeRecordConstructor::records->at(RootRecord::record_offset+ROOT_RECORD_SIZE+index)
        = offset;

}


size_t SnarlDistanceIndex::SnarlRecord::distance_vector_size(record_t type, size_t node_count) {
    if (type == SNARL || type == ROOT_SNARL){
        //For a normal snarl, its just the record size and the pointers to children
        return 0;
    } else if (type == DISTANCED_SNARL) {
        //For a normal min distance snarl, record size and the pointers to children, and
        //matrix of distances
        size_t node_side_count = node_count * 2 + 2;
        return  (((node_side_count+1)*node_side_count) / 2);
    } else if (type ==  OVERSIZED_SNARL){
        //For a large min_distance snarl, record the side, pointers to children, and just
        //the min distances from each node side to the two boundary nodes
        size_t node_side_count = node_count * 2 + 2;
        return  (node_side_count * 2);
    } else if (type == DISTANCED_ROOT_SNARL) {
        size_t node_side_count = node_count * 2;
        return  (((node_side_count+1)*node_side_count) / 2);
    } else {
        throw runtime_error ("error: this is not a snarl");
    }
}

size_t SnarlDistanceIndex::SnarlRecord::record_size (record_t type, size_t node_count) {
    return SNARL_RECORD_SIZE + distance_vector_size(type, node_count) + node_count;
}
size_t SnarlDistanceIndex::SnarlRecord::record_size() {
    record_t type = get_record_type();
    return record_size(type, get_node_count());
}


int64_t SnarlDistanceIndex::SnarlRecord::get_distance_vector_offset(size_t rank1, bool right_side1, size_t rank2,
        bool right_side2, size_t node_count, record_t type) {

    //how many node sides in this snarl
    size_t node_side_count = type == DISTANCED_ROOT_SNARL ? node_count * 2 : node_count * 2 + 2;

    //make sure we're looking at the correct node side
    //If this is the start or end node, then we don't adjust based on orientation
    //because we only care about the inner sides. If this is a root snarl, then
    //there is no start or end node and the ranks 0 and 1 are not special
     if (type == DISTANCED_ROOT_SNARL) {
        rank1 = rank1 * 2;
        if (right_side1) {
            rank1 += 1;
        }
     } else if (rank1 != 0 && rank1 != 1) {
        rank1 = (rank1-1) * 2;
        if (right_side1) {
            rank1 += 1;
        }
    }
    if (type == DISTANCED_ROOT_SNARL) {
        rank2 = rank2 * 2;
        if (right_side2) {
            rank2 += 1;
        }
    } else if (rank2 != 0 && rank2 != 1) {
        rank2 = (rank2-1) * 2;
        if (right_side2) {
            rank2 += 1;
        }
    }

    //reverse order of ranks if necessary
    if (rank1 > rank2) {
        size_t tmp = rank1;
        rank1 = rank2;
        rank2 = tmp;
    }

    if (type == SNARL || type == ROOT_SNARL) {
        throw runtime_error("error: trying to access distance in a distanceless snarl tree");
    } else if (type == DISTANCED_SNARL || type == DISTANCED_ROOT_SNARL) {
        //normal distance index
        size_t k = node_side_count-rank1;
        return (((node_side_count+1) * node_side_count)/2) - (((k+1)*k) / 2) + rank2 - rank1;
    } else if (type ==  OVERSIZED_SNARL) {
        //abbreviated distance index storing only the distance from each node side to the
        //start and end
        if (rank1 == 0) {
            return rank2;
        } else if (rank1 == 1) {
            return node_count + 2 + rank2;
        } else {
            throw runtime_error("error: trying to access distance in an oversized snarl");
        }
    } else {
        throw runtime_error("error: trying to distance from something that isn't a snarl");
    }
}

int64_t SnarlDistanceIndex::SnarlRecord::get_distance_vector_offset(size_t rank1, bool right_side1,
        size_t rank2, bool right_side2) const {
    return get_distance_vector_offset(rank1, right_side1, rank2, right_side2,
        get_node_count(), get_record_type());
}

int64_t SnarlDistanceIndex::SnarlRecord::get_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2) const {

    //offset of the start of the distance vectors in snarl_tree_records
    size_t distance_vector_start = record_offset + SNARL_RECORD_SIZE + get_node_count();
    //Offset of this particular distance in the distance vector
    size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

    size_t val = records->at(distance_vector_start + distance_vector_offset);
    return  val == 0 ? std::numeric_limits<int64_t>::max() : val-1;

}

bool SnarlDistanceIndex::SnarlRecord::for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {
    size_t child_count = get_node_count();
    size_t child_record_offset = get_child_record_pointer();
    for (size_t i = 0 ; i < child_count ; i++) {
        size_t child_offset =  records->at(child_record_offset + i);
        net_handle_record_t type = SnarlTreeRecord(child_offset, records).get_record_handle_type();
        assert(type == NODE_HANDLE || type == CHAIN_HANDLE);
        net_handle_t child_handle =  get_net_handle (child_offset, START_END, CHAIN_HANDLE);
        bool result = iteratee(child_handle);
        if (result == false) {
            return false;
        }
    }
    return true;
}


void SnarlDistanceIndex::SnarlRecordConstructor::set_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2, int64_t distance) {
    //offset of the start of the distance vectors in snarl_tree_records
    size_t distance_vector_start = SnarlRecord::record_offset + SNARL_RECORD_SIZE + get_node_count();
    //Offset of this particular distance in the distance vector
    size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

    size_t val = distance == std::numeric_limits<int64_t>::max() ? 0 : distance+1;
#ifdef debug_indexing
    cerr <<  distance_vector_start + distance_vector_offset << " set distance_value " << val << endl;
    //assert(SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) == 0 ||
    //        SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) == val);
#endif


    SnarlTreeRecordConstructor::records->at(distance_vector_start + distance_vector_offset) = val;
}

//Node count is the number of nodes in the snarl, not including boundary nodes
void SnarlDistanceIndex::SnarlRecordConstructor::set_node_count(size_t node_count) {
#ifdef debug_indexing
    cerr << SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET << " set snarl node count " << node_count << endl;
    assert(node_count > 0);//TODO: Don't bother making a record for trivial snarls
    assert(SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset + SNARL_NODE_COUNT_OFFSET) = node_count;
}


void SnarlDistanceIndex::SnarlRecordConstructor::set_child_record_pointer(size_t pointer) {
#ifdef debug_indexing
    cerr << SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET << " set snarl child record pointer " << pointer << endl;
    assert(SnarlTreeRecordConstructor::SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(SnarlRecord::record_offset+SNARL_CHILD_RECORD_OFFSET) = pointer;
}
//Add a reference to a child of this snarl. Assumes that the index is completed up
//to here
void SnarlDistanceIndex::SnarlRecordConstructor::add_child(size_t pointer){
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << " Adding child pointer to the end of the array " << endl;
#endif
    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+1);
    SnarlTreeRecordConstructor::records->at(start_i) = pointer;
}

size_t SnarlDistanceIndex::ChainRecord::get_node_count() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        return 1;
    } else {
        return records->at(record_offset + CHAIN_NODE_COUNT_OFFSET);
    }
}

pair<size_t, bool> SnarlDistanceIndex::ChainRecord::get_last_child_offset() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get children of a node");
    } else {
        size_t val = records->at(record_offset + CHAIN_LAST_CHILD_OFFSET) < 1;
        return make_pair(val<1, val & 1);
    }
}
bool SnarlDistanceIndex::ChainRecord::get_is_looping_chain_connected_backwards() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get children of a node");
    } else {
        size_t val = records->at(record_offset + CHAIN_LAST_CHILD_OFFSET);
        return val & 1;
    }
}

int64_t SnarlDistanceIndex::ChainRecord::get_chain_node_id(size_t pointer) const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    }
    return records->at(pointer+CHAIN_NODE_ID_OFFSET);
}
int64_t SnarlDistanceIndex::ChainRecord::get_prefix_sum_value(size_t pointer) const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    }
    size_t val = records->at(pointer+CHAIN_NODE_PREFIX_SUM_OFFSET);
    return val == 0 ? std::numeric_limits<int64_t>::max() : val-1;
}
int64_t SnarlDistanceIndex::ChainRecord::get_forward_loop_value(size_t pointer) const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    }
    size_t val = records->at(pointer+CHAIN_NODE_FORWARD_LOOP_OFFSET);
    return val == 0 ? std::numeric_limits<int64_t>::max() : val-1;
}
int64_t SnarlDistanceIndex::ChainRecord::get_reverse_loop_value(size_t pointer) const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    }
    size_t val =  records->at(pointer+CHAIN_NODE_REVERSE_LOOP_OFFSET);
    return val == 0 ? std::numeric_limits<int64_t>::max() : val-1;
}
size_t SnarlDistanceIndex::ChainRecord::get_chain_component(size_t pointer, bool get_end) const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    }
    if (get_record_type() != MULTICOMPONENT_CHAIN){
        throw runtime_error("error: Trying to get the component of a single component chain");
    }
    if (get_chain_node_id(pointer) == get_start_id() && !get_end) {
        //The first component is always 0
        return 0;
    } else {
        return records->at(pointer+CHAIN_NODE_COMPONENT_OFFSET);
    }
}

int64_t SnarlDistanceIndex::ChainRecord::get_distance(tuple<size_t, bool, int64_t> node1,
                     tuple<size_t, bool, int64_t> node2) const {

    if (std::get<0>(node1) > std::get<0>(node2)) {
        //If the first node comes after the second in the chain, reverse them
        tuple<size_t, bool, size_t> tmp = node1;
        node1 = node2;
        node2 = tmp;

    }

    bool is_looping_chain = get_start_id() == get_end_id();
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    } else if (get_record_type() == MULTICOMPONENT_CHAIN) {
        if (get_chain_component(std::get<0>(node1)) != get_chain_component(std::get<0>(node2))) {
            if (is_looping_chain) {
                //If this is a looping chain, then the first/last node could be in two
                //components
                return get_distance_taking_chain_loop(node1, node2);
            } else {
                return std::numeric_limits<int64_t>::max();
            }
        }
    }


    int64_t distance;

    if (!std::get<1>(node1) && std::get<1>(node2)) {
        //Right of 1 and left of 2, so a simple forward traversal of the chain
        if (std::get<0>(node1) == std::get<0>(node2)) {
            //If these are the same node, then the path would need to go around the node
            distance = sum({get_forward_loop_value(std::get<0>(node1)),
                        get_reverse_loop_value(std::get<0>(node1)),
                        std::get<2>(node1)});
        } else {
            distance = minus(get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)),
                 std::get<2>(node1));
        }
    } else if (!std::get<1>(node1) && !std::get<1>(node2)) {
        //Right side of 1 and right side of 2
        if (std::get<0>(node1) == std::get<0>(node2)) {
            distance = get_forward_loop_value(std::get<0>(node2));

        } else {
            distance = minus( sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)) ,
                               std::get<2>(node2),
                               get_forward_loop_value(std::get<0>(node2))}),
                         std::get<2>(node1));
        }
    } else if (std::get<1>(node1) && std::get<1>(node2)) {
        //Left side of 1 and left side of 2
        if (std::get<0>(node1) == std::get<0>(node2)) {
            distance = get_reverse_loop_value(std::get<0>(node1));

        } else {
            distance = sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)),
                            get_reverse_loop_value(std::get<0>(node1))});
        }
    } else {
        assert(std::get<1>(node1) && !std::get<1>(node2));
        //Left side of 1 and right side of 2
        distance = sum({get_prefix_sum_value(std::get<0>(node2)) - get_prefix_sum_value(std::get<0>(node1)),
                        get_reverse_loop_value(std::get<0>(node1)),
                        get_forward_loop_value(std::get<0>(node2)),
                        std::get<2>(node2)});

    }
    if (is_looping_chain) {
        distance = std::min(distance, get_distance_taking_chain_loop(node1, node2));
    }
    return distance;
}


int64_t SnarlDistanceIndex::ChainRecord::get_distance_taking_chain_loop(tuple<size_t, bool, int64_t> node1,
                     tuple<size_t, bool, int64_t> node2) const {
    //This is only called by get_distance, so the nodes should be ordered
    assert (std::get<0>(node1) <= std::get<0>(node2));

    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain distances from a node");
    } else if (get_record_type() == MULTICOMPONENT_CHAIN) {
        size_t last_component = get_chain_component(get_first_node_offset(), true);
        bool first_in_first_component = get_chain_component(std::get<0>(node1)) == 0 || get_chain_component(std::get<0>(node1)) == last_component;
        bool second_in_first_component = get_chain_component(std::get<0>(node2)) == 0 || get_chain_component(std::get<0>(node2)) == last_component;
        bool can_loop = get_is_looping_chain_connected_backwards();

        if (!can_loop || !first_in_first_component || ! second_in_first_component) {
            //If this is a multicomponent chain then it can only reach around backwards if both nodes
            //are in the first (last) component
            return std::numeric_limits<int64_t>::max();
        }
    }


    int64_t distance;
    assert(get_start_id() == get_end_id());

    if (!std::get<1>(node1) && std::get<1>(node2)) {
        //Right of 1 and left of 2, so a simple forward traversal of the chain
        //loop forward from the first node, from the start of the chain to the first
        //node, from the end of the node to the second node, and the reverse loop of the second
        distance = sum({get_forward_loop_value(std::get<0>(node1)),
                             std::get<2>(node1),
                             get_prefix_sum_value(std::get<0>(node1)),
                             minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                             get_reverse_loop_value(std::get<0>(node2))});
    } else if (!std::get<1>(node1) && !std::get<1>(node2)) {
        //Right side of 1 and right side of 2

        //Check distance for taking loop in chain: loop forward from the first node, from the start of the
        //chain to the first node, from the end of the node to the second node
        distance = sum({get_forward_loop_value(std::get<0>(node1)),
                             std::get<2>(node1),
                             get_prefix_sum_value(std::get<0>(node1)),
                             minus(minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                                std::get<2>(node2))});
    } else if (std::get<1>(node1) && std::get<1>(node2)) {
        //Left side of 1 and left side of 2

        //from the first node left to the start, around the
        //chain loop, then the reverse loop of the second node
        //This assumes that the length of the chain only includes the start/end node's length once,
        //which it does but might change
        distance = sum({get_prefix_sum_value(std::get<0>(node1)),
                             minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                             get_reverse_loop_value(std::get<0>(node2))});
    } else {
        assert(std::get<1>(node1) && !std::get<1>(node2));
        //Left side of 1 and right side of 2

        //Check the distance going backwards around the chain
        distance = sum({get_prefix_sum_value(std::get<0>(node1)),
                             minus(minus(get_min_length(), get_prefix_sum_value(std::get<0>(node2))),
                              std::get<2>(node2))});
    }
    return distance;
}

size_t SnarlDistanceIndex::ChainRecord::get_first_node_offset() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain traversal from a node");
    }
    return record_offset + CHAIN_RECORD_SIZE;
}

pair<size_t, bool> SnarlDistanceIndex::ChainRecord::get_next_child(const pair<size_t, bool> pointer, bool go_left) const {
    //If this is a multicomponent chain, then the size of each node record in the chain is bigger
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain traversal from a node");
    }
    size_t node_record_size = get_record_type() == MULTICOMPONENT_CHAIN ? CHAIN_NODE_MULTICOMPONENT_RECORD_SIZE : CHAIN_NODE_RECORD_SIZE;
    if (pointer.second) {
        //This is a snarl
        if (go_left) {
            return make_pair(pointer.first - node_record_size, false);
        } else {
            if (get_start_id() == get_end_id() && pointer.first == get_last_child_offset().first) {
                //If this is the last child in a looping chain
                return make_pair(get_first_node_offset(), false);
            } else {
                size_t snarl_record_length = SnarlRecord(pointer.first, records).record_size();
                return make_pair(pointer.first + snarl_record_length + 1, false);
            }
        }
    } else {
        //This is a node
        if (go_left) {
            //Looking left in the chain
            if (pointer.first == get_first_node_offset()) {
                //If this is the first node in the chain
                if (get_start_id() == get_end_id()) {
                    pair<size_t, bool> last_child = get_last_child_offset();
                    return last_child; //TODO: I'm not sure if this is always true (snarl)
                } else {
                    return make_pair(0, false);
                }
            }
            size_t snarl_record_size = records->at(pointer.first-1);
            if (snarl_record_size == 0) {
                //Just another node to the left
                return make_pair(pointer.first-node_record_size, false);
            } else {
                //There is a snarl to the left of this node
                return make_pair(pointer.first - snarl_record_size - 1, true);
            }
        } else {
            //Looking right in the chain
            if (records->at(pointer.first+node_record_size-1) == 0 &&
                records->at(pointer.first+node_record_size) == 0) {
                //If this is the last node in the chain
                if (get_start_id() == get_end_id()) {
                    //If this loops, go back to the beginning
                    return make_pair(get_first_node_offset(), false);
                } else {
                    return make_pair(0, false);
                }
            }
            size_t snarl_record_size = records->at(pointer.first+node_record_size-1);
            return make_pair(pointer.first+node_record_size, snarl_record_size != 0);
        }
    }
}
net_handle_t SnarlDistanceIndex::ChainRecord::get_next_child(const net_handle_t& net_handle, bool go_left) const {
    //get the next child in the chain. net_handle must point to a snarl or node in the chain
    net_handle_record_t handle_type = get_handle_type(net_handle);
    net_handle_record_t record_type = get_record_handle_type();
    if (record_type == NODE_HANDLE) {
        //If this record is a node pretending to be a chain, then there is no next child
        assert(handle_type == NODE_HANDLE);
        assert(get_record_offset(net_handle) == record_offset);
        return net_handle;
    }

    //Get the current pointer, pointing at the net_handle in the chain
    if (handle_type == NODE_HANDLE) {
        //If this net handle is a node, then it's rank in parent points to it in the chain
        NodeRecord node_record(get_record_offset(net_handle), records);
        assert(node_record.get_parent_record_offset() == record_offset);
        pair<size_t, bool> next_pointer = get_next_child(
            make_pair(node_record.get_rank_in_parent(), false), go_left);
        if (next_pointer.first == 0 ){
            return net_handle;
        }
        bool next_is_reversed_in_parent = NodeRecord(
                get_offset_from_id(records->at(next_pointer.first)), records
            ).get_is_reversed_in_parent();

        connectivity_t connectivity = go_left == next_is_reversed_in_parent ? START_END : END_START;
        if (!next_pointer.second && next_pointer.first == get_first_node_offset()) {
            connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), get_end_endpoint(connectivity));
        }
        return get_net_handle(get_offset_from_id(records->at(next_pointer.first)),
                          connectivity,
                          next_pointer.second ? SNARL_HANDLE : NODE_HANDLE);
    } else {
        //Otherwise, it is a snarl and we can use the snarl's offset, since it exists in
        //the chain
        assert(handle_type == SNARL_HANDLE) ;
        pair<size_t, bool> next_pointer = get_next_child(
            make_pair(get_record_offset(net_handle), true), go_left);
        if (next_pointer.first == 0 ){
            return net_handle;
        }
        connectivity_t connectivity = next_pointer.second ? END_START : START_END;
        if (!next_pointer.second && next_pointer.first == get_first_node_offset()) {
            connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), get_end_endpoint(connectivity));
        }
        return get_net_handle(next_pointer.first, connectivity,
                          next_pointer.first ? SNARL_HANDLE : NODE_HANDLE);
    }
}

bool SnarlDistanceIndex::ChainRecord::for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {

    if (get_record_handle_type() == NODE_HANDLE) {
        //If this is a node pretending to be a chain, just do it for the node
        return iteratee(get_net_handle(record_offset, START_END, NODE_HANDLE));
    }


    //If this is a node, then the offset of the node in the chain, false
    //If it is a snarl, then the offset of the snarl record, true
    pair<size_t, bool> current_child (get_first_node_offset(), false);
    bool is_first = true;

    while (current_child.first != 0) {
        net_handle_t child_handle = current_child.second 
            ? get_net_handle (current_child.first, START_END, SNARL_HANDLE)
            : get_net_handle (get_offset_from_id(records->at(current_child.first)), START_END, NODE_HANDLE);
        if (!is_first && current_child == make_pair(get_first_node_offset(), false)){
            //Don't look at the first node a second time
            return true;
        }

        bool result = iteratee(child_handle);
        if (result == false) {
            return false;
        }
        current_child = get_next_child(current_child, false);
        is_first = false;
    }
    return true;
}


void SnarlDistanceIndex::ChainRecordConstructor::add_node(id_t id, int64_t prefix_sum, int64_t forward_loop, int64_t reverse_loop) {
    assert(ChainRecord::get_record_type() != MULTICOMPONENT_CHAIN);
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << " - " << SnarlTreeRecordConstructor::records->size() + 3 << " Adding chain's child node " << id << " to the end of the a     rray (values: " << prefix_sum << "," << forward_loop << "," << reverse_loop << ")" << endl;
#endif

    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+4);

    SnarlTreeRecordConstructor::records->at(start_i) = id;
    SnarlTreeRecordConstructor::records->at(start_i+1) = 
            prefix_sum==std::numeric_limits<int64_t>::max() ? 0 : prefix_sum+1;
    SnarlTreeRecordConstructor::records->at(start_i+2) = 
            forward_loop==std::numeric_limits<int64_t>::max() ? 0 : forward_loop+1;
    SnarlTreeRecordConstructor::records->at(start_i+3) = 
            reverse_loop==std::numeric_limits<int64_t>::max() ? 0 : reverse_loop+1;
}
void SnarlDistanceIndex::ChainRecordConstructor::add_node(id_t id, int64_t prefix_sum, int64_t forward_loop, int64_t reverse_loop, size_t component) {
    assert(ChainRecord::get_record_type() == MULTICOMPONENT_CHAIN);
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << " - " << SnarlTreeRecordConstructor::records->size() + 4 << " Adding chain's child node the end of the array " << endl;
#endif

    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+5);

    SnarlTreeRecordConstructor::records->at(start_i) = id;
    SnarlTreeRecordConstructor::records->at(start_i+1) = prefix_sum==std::numeric_limits<int64_t>::max() ? 0 : prefix_sum+1;
    SnarlTreeRecordConstructor::records->at(start_i+2) = forward_loop==std::numeric_limits<int64_t>::max() ? 0 : forward_loop+1;
    SnarlTreeRecordConstructor::records->at(start_i+3) = reverse_loop==std::numeric_limits<int64_t>::max() ? 0 : reverse_loop+1;
    SnarlTreeRecordConstructor::records->at(start_i+4) = component;
}
void SnarlDistanceIndex::ChainRecordConstructor::set_node_count(size_t node_count) {
#ifdef debug_indexing
    cerr << ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET << " set chain node count " << node_count << endl;
    assert(SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_NODE_COUNT_OFFSET) = node_count;
}

//The offset of the last child, if it is a snarl, and if it can loop
void SnarlDistanceIndex::ChainRecordConstructor::set_last_child_offset(size_t offset, bool is_snarl, bool loopable) {
#ifdef debug_indexing
    cerr << ChainRecord::record_offset + CHAIN_LAST_CHILD_OFFSET << " set chain last child offset " << offset << endl;
    assert(SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_LAST_CHILD_OFFSET) == 0);
#endif

    SnarlTreeRecordConstructor::records->at(ChainRecord::record_offset + CHAIN_LAST_CHILD_OFFSET) = ((offset < 2) | (is_snarl<1)) | loopable;
}

void SnarlDistanceIndex::ChainRecordConstructor::add_trivial_snarl() {
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << "  Adding chain's trivial snarl to the end of the array " << endl;
#endif
    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+1);

    SnarlTreeRecordConstructor::records->at(start_i) = 0;
}
//Add a snarl to the end of the chain and return a SnarlRecordConstructor pointing to it
SnarlDistanceIndex::SnarlRecordConstructor SnarlDistanceIndex::ChainRecordConstructor::add_snarl(size_t snarl_size, record_t type) {
    size_t snarl_record_size = SnarlRecord::record_size(type, snarl_size);
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << " Adding child snarl length to the end of the array " << endl;
#endif

    
    
    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+1);
    SnarlTreeRecordConstructor::records->at(start_i) = snarl_record_size;
    SnarlTreeRecordConstructor::records->reserve(start_i + snarl_record_size);
    SnarlRecordConstructor snarl_record(snarl_size, SnarlTreeRecordConstructor::records, type);
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size() << " Adding child snarl length to the end of the array " << endl;
#endif
    start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+1);
    SnarlTreeRecordConstructor::records->at(start_i) = snarl_record_size;
    return snarl_record;
}
void SnarlDistanceIndex::ChainRecordConstructor::finish_chain(){
#ifdef debug_indexing
    cerr << SnarlTreeRecordConstructor::records->size()  << " - " <<  SnarlTreeRecordConstructor::records->size()+1 << " Adding the last two chain 0's to the end of the array " <<      endl;  
#endif

    size_t start_i = SnarlTreeRecordConstructor::records->size();
    SnarlTreeRecordConstructor::records->resize(start_i+2);
    SnarlTreeRecordConstructor::records->at(start_i) = 0;
    SnarlTreeRecordConstructor::records->at(start_i+1) = 0;
}
string SnarlDistanceIndex::net_handle_as_string(const net_handle_t& net) const {
    net_handle_record_t type = get_handle_type(net);
    SnarlTreeRecord record (net, &snarl_tree_records);
    net_handle_record_t record_type = record.get_record_handle_type();
    string result;
    if (type == ROOT_HANDLE) {
        return "root"; 
    } else if (type == NODE_HANDLE) {
        if (ends_at(net) == starts_at(net)) {
            return "node" + std::to_string( get_node_id_from_offset(get_record_offset(net))) + (ends_at(net) == START ? "rev" : "fd") + " that is the end node of a looping chain";
        }
        return  "node " + std::to_string( get_node_id_from_offset(get_record_offset(net))) + (ends_at(net) == START ? "rev" : "fd");
    } else if (type == SNARL_HANDLE) {
        if (record.get_record_type() == ROOT) {
            return "root snarl";
        }
        result += "snarl ";         
    } else if (type == CHAIN_HANDLE && record_type == NODE_HANDLE) {
        return  "node " + std::to_string( get_node_id_from_offset(get_record_offset(net)))
               + (ends_at(net) == START ? "rev" : "fd") + " pretending to be a chain";
    } else if (type == CHAIN_HANDLE) {
        result += "chain ";
    } else if (type == SENTINEL_HANDLE) {
        result += "sentinel of snarl ";
    } else {
        throw runtime_error("error: Unknown net_handle_t type");
    }
    result += ( std::to_string(record.get_start_id())
            + (record.get_start_orientation() ? "rev" : "fd")
            + "->"
            + std::to_string(record.get_end_id())
            + (record.get_end_orientation() ? "rev" : "fd"));
    result += "traversing ";
    result += (starts_at(net) == START ? "start" : (starts_at(net) == END ? "end" : "tip"));
    result += "->";
    result += (ends_at(net) == START ? "start" : (ends_at(net) == END ? "end" : "tip"));
    return result;
}


}
