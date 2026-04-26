//#define debug_distance_indexing
//#define debug_snarl_traversal
//#define debug_distances
//#define debug_hub_label_build
//#define debug_hub_label_storage

#include "snarl_distance_index.hpp"
#include "snarl_distance_index_child_graph.hpp"
#include <span>

using namespace std;
using namespace handlegraph;
namespace vg {

SnarlDistanceIndex::TemporaryDistanceIndex make_temporary_distance_index(
    const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit, bool only_top_level_chain_distances)  {

#ifdef debug_distance_indexing
    cerr << "Creating new distance index for nodes between " << graph->min_node_id() << " and " << graph->max_node_id() << endl;

#endif

    SnarlDistanceIndex::TemporaryDistanceIndex temp_index;

    temp_index.min_node_id=graph->min_node_id();
    temp_index.max_node_id=graph->max_node_id();

    //Construct the distance index using the snarl decomposition
    //traverse_decomposition will visit all structures (including trivial snarls), calling
    //each of the given functions for the start and ends of the snarls and chains

    temp_index.temp_node_records.resize(temp_index.max_node_id-temp_index.min_node_id+1);



    //Stores unfinished records, as type of record and offset into appropriate vector
    //(temp_node/snarl/chain_records)
    vector<SnarlDistanceIndex::temp_record_ref_t> stack;

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
        //assert(temp_index.get_node(make_pair(SnarlDistanceIndex::TEMP_NODE, graph->get_id(chain_start_handle))).node_id == 0);
#endif

        //Fill in node in chain
        stack.emplace_back(SnarlDistanceIndex::TEMP_CHAIN, temp_index.temp_chain_records.size());
        nid_t node_id = graph->get_id(chain_start_handle);
        temp_index.temp_chain_records.emplace_back();
        auto& temp_chain = temp_index.temp_chain_records.back();
        temp_chain.start_node_id = node_id; 
        temp_chain.start_node_rev = graph->get_is_reverse(chain_start_handle);
        temp_chain.children.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);


        //And the node record itself
        auto& temp_node = temp_index.get_node(temp_chain.children.back());
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
        SnarlDistanceIndex::temp_record_ref_t chain_index = stack.back();
        stack.pop_back();

#ifdef debug_distance_indexing
        assert(chain_index.first == SnarlDistanceIndex::TEMP_CHAIN);
#endif
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = temp_index.get_chain(chain_index);
        nid_t node_id = graph->get_id(chain_end_handle);

        if (temp_chain_record.children.size() == 1 && node_id == temp_chain_record.start_node_id) {
            //This is a trivial snarl

#ifdef debug_distance_indexing
            //Then this must be the last thing on the chain_records vector
            assert(temp_index.temp_chain_records.size() == chain_index.second+1);
#endif

            //Get the node
            SnarlDistanceIndex::temp_record_ref_t node_index = make_pair(SnarlDistanceIndex::TEMP_NODE, node_id);
            SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index.get_node(node_index);

            temp_node_record.reversed_in_parent = false;

            //And give the chain's parent the node info
            //
            if (stack.empty()) {
                temp_node_record.parent = make_pair(SnarlDistanceIndex::TEMP_ROOT, 0);
                //If this was the last thing on the stack, then this was a root

                //Check to see if there is anything connected to the ends of the chain
                vector<nid_t> reachable_nodes;
                graph->follow_edges(graph->get_handle(node_id, false),
                    false, [&] (const handle_t& next) {
                        if (graph->get_id(next) != node_id) {
                            reachable_nodes.emplace_back(graph->get_id(next));
                        }
                    });
                graph->follow_edges(graph->get_handle(node_id, true),
                    false, [&] (const handle_t& next) {
                        if (graph->get_id(next) != node_id) {
                            reachable_nodes.emplace_back(graph->get_id(next));
                        }
                    });
                if (reachable_nodes.size()) {
                    //If we can reach anything leaving the chain (besides the chain itself), then it is part of a root snarl
                    //Note that if the chain's start and end node are the same, then it will always be a single component
#ifdef debug_distance_indexing
                    cerr << "                 This trivial chain is part of the root but connects with something else in the root"<<endl;
#endif
                    bool new_component = true;

                    //Add this to the union find
                    root_snarl_component_uf.resize(root_snarl_component_uf.size() + 1);
                    //And remember that it's in a connected component of the root
                    temp_node_record.root_snarl_index = temp_index.root_snarl_components.size();
                    temp_index.root_snarl_components.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);
                    for (nid_t next_id : reachable_nodes) {
                        //For each node that this is connected to, check if we've already seen it and if we have, then
                        //union this chain and that node's chain
                        SnarlDistanceIndex::temp_record_ref_t next_index = make_pair(SnarlDistanceIndex::TEMP_NODE, next_id);
                        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& node_record = temp_index.get_node(next_index);
                        if (node_record.node_id != 0) {
                            //If we've already seen this node, union it with the new one
                            //If we can see it by walking out from this top-level chain, then it must also be a
                            //top-level chain (or node pretending to be a chain)
                            size_t other_i = node_record.parent.first == SnarlDistanceIndex::TEMP_CHAIN
                                           ? temp_index.get_chain(node_record.parent).root_snarl_index
                                           : node_record.root_snarl_index;
#ifdef debug_distance_indexing
                            assert(other_i != std::numeric_limits<size_t>::max());
#endif
                            root_snarl_component_uf.union_groups(other_i, temp_node_record.root_snarl_index);
//#ifdef debug_distance_indexing
//                            cerr << "        Union this trivial  with " << temp_index.get_chain(node_record.parent).start_node_id << " " << temp_index.get_chain(node_record.parent).end_node_id << endl;
//#endif
                        } else {
                            new_component = false;
                        }
                    }
                } else {
                    //If this chain isn't connected to anything else, then it is a single component of the root
                    temp_node_record.rank_in_parent = temp_index.components.size();
                    temp_index.components.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);
                }
            } else {
                //The last thing on the stack is the parent of this chain, which must be a snarl
                temp_node_record.parent = stack.back();
                auto& parent_snarl_record = temp_index.get_snarl(temp_node_record.parent);
                temp_node_record.rank_in_parent = parent_snarl_record.children.size() + 2;
                parent_snarl_record.children.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);
            }


            //Remove the chain record
            temp_index.temp_chain_records.pop_back();
            temp_index.max_index_size += temp_node_record.get_max_record_length();

        } else {
            //Otherwise, it is an actual chain

            //Fill in node in chain
            temp_chain_record.end_node_id = node_id;
            temp_chain_record.end_node_rev = graph->get_is_reverse(chain_end_handle);
            temp_chain_record.end_node_length = graph->get_length(chain_end_handle);
            
            bool is_root_chain = false;

            if (stack.empty()) {
                //If this was the last thing on the stack, then this was a root
                is_root_chain = true;

                //Check to see if there is anything connected to the ends of the chain
                vector<nid_t> reachable_nodes;
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
                    temp_chain_record.root_snarl_index = temp_index.root_snarl_components.size();
                    temp_index.root_snarl_components.emplace_back(chain_index);
                    for (nid_t next_id : reachable_nodes) {
                        //For each node that this is connected to, check if we've already seen it and if we have, then
                        //union this chain and that node's chain
                        SnarlDistanceIndex::temp_record_ref_t next_index = make_pair(SnarlDistanceIndex::TEMP_NODE, next_id);
                        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& node_record = temp_index.get_node(next_index);
                        if (node_record.node_id != 0) {
                            //If we've already seen this node, union it with the new one
                            //If we can see it by walking out from this top-level chain, then it must also be a
                            //top-level chain (or node pretending to be a chain)
                            size_t other_i = node_record.parent.first == SnarlDistanceIndex::TEMP_CHAIN
                                           ? temp_index.get_chain(node_record.parent).root_snarl_index
                                           : node_record.root_snarl_index;
#ifdef debug_distance_indexing
                            assert(other_i != std::numeric_limits<size_t>::max());
#endif
                            root_snarl_component_uf.union_groups(other_i, temp_chain_record.root_snarl_index);
#ifdef debug_distance_indexing
                            if (node_record.parent.first == SnarlDistanceIndex::TEMP_CHAIN) {
                                cerr << "        Union this chain with " << temp_index.get_chain(node_record.parent).start_node_id << " " << temp_index.get_chain(node_record.parent).end_node_id << endl;
                            } else {
                                cerr << "        Union this chain with root " << node_record.root_snarl_index << endl;
                            }
#endif
                        } else {
                            new_component = false;
                        }
                    }
                } else {
                    //If this chain isn't connected to anything else, then it is a single component of the root
                    temp_chain_record.parent = make_pair(SnarlDistanceIndex::TEMP_ROOT, 0);
                    temp_chain_record.rank_in_parent = temp_index.components.size();
                    temp_index.components.emplace_back(chain_index);
                }
            } else {
                //The last thing on the stack is the parent of this chain, which must be a snarl
                temp_chain_record.parent = stack.back();
                auto& parent_snarl_record = temp_index.get_snarl(temp_chain_record.parent);
                temp_chain_record.rank_in_parent = parent_snarl_record.children.size() + 2;
                parent_snarl_record.children.emplace_back(chain_index);
            }

            temp_index.max_index_size += temp_chain_record.get_max_record_length(!only_top_level_chain_distances || is_root_chain ? true : false );
#ifdef debug_distance_indexing
            cerr << "  Ending new " << (temp_chain_record.is_trivial ? "trivial " : "") <<  "chain " << temp_index.structure_start_end_as_string(chain_index)
              << endl << "    that is a child of " << temp_index.structure_start_end_as_string(temp_chain_record.parent) << endl;
#endif
        }
    },
    [&](handle_t snarl_start_handle) {
        /*This gets called at the beginning of a new snarl facing in
         * Create a new snarl record and fill in the start node.
         * The node record would have been created as part of the chain, or as the end node
         * of the previous snarl
         */

#ifdef debug_distance_indexing
        cerr << "  Starting new snarl at " << graph->get_id(snarl_start_handle) << (graph->get_is_reverse(snarl_start_handle) ? " reverse" : " forward") << endl;
        cerr << "with index " << temp_index.temp_snarl_records.size() << endl;
#endif
        auto& parent = stack.back();
        stack.emplace_back(SnarlDistanceIndex::TEMP_SNARL, temp_index.temp_snarl_records.size());
        temp_index.temp_snarl_records.emplace_back();
        temp_index.temp_snarl_records.back().start_node_id = graph->get_id(snarl_start_handle);
        temp_index.temp_snarl_records.back().start_node_rev = graph->get_is_reverse(snarl_start_handle);
        temp_index.temp_snarl_records.back().start_node_length = graph->get_length(snarl_start_handle);

    },
    [&](handle_t snarl_end_handle){
        /*This gets called at the end of the snarl facing out
         * Fill in the end node of the snarl, its parent, and record the snarl as a child of its
         * parent chain
         * Also create a node record
         */
        SnarlDistanceIndex::temp_record_ref_t snarl_index = stack.back();
        stack.pop_back();
#ifdef debug_distance_indexing
        assert(snarl_index.first == SnarlDistanceIndex::TEMP_SNARL);
        assert(stack.back().first == SnarlDistanceIndex::TEMP_CHAIN);
#endif
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(snarl_index);
        nid_t node_id = graph->get_id(snarl_end_handle);

        //Record the end node in the snarl
        temp_snarl_record.end_node_id = node_id;
        temp_snarl_record.end_node_rev = graph->get_is_reverse(snarl_end_handle);
        temp_snarl_record.end_node_length = graph->get_length(snarl_end_handle);
        temp_snarl_record.node_count = temp_snarl_record.children.size();
        bool any_edges_in_snarl = false;
        graph->follow_edges(graph->get_handle(temp_snarl_record.start_node_id, temp_snarl_record.start_node_rev), false, [&](const handle_t& next_handle) {
            if (graph->get_id(next_handle) != temp_snarl_record.end_node_id) {
                any_edges_in_snarl = true;
            }
        });
        graph->follow_edges(graph->get_handle(temp_snarl_record.end_node_id, !temp_snarl_record.end_node_rev), false, [&](const handle_t& next_handle) {
            if (graph->get_id(next_handle) != temp_snarl_record.start_node_id) {
                any_edges_in_snarl = true;
            }
        });

        if (temp_snarl_record.children.size() == 0) {
            //This is a trivial snarl
            temp_snarl_record.is_trivial = true;

#ifdef debug_distance_indexing
            cerr << "  Ending and forgetting trivial snarl " << temp_index.structure_start_end_as_string(snarl_index)
                 << endl << "    that is a child of " << temp_index.structure_start_end_as_string(temp_snarl_record.parent) << endl;
#endif

            //Add the end node to the chain
#ifdef debug_distance_indexing
            assert(stack.back().first == SnarlDistanceIndex::TEMP_CHAIN);
#endif
            temp_snarl_record.parent = stack.back();
            auto& temp_chain = temp_index.get_chain(stack.back());
            temp_chain.children.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);

            //Remove the snarl record.
            //This invalidates snarl_index!!!
#ifdef debug_distance_indexing
            assert(temp_index.temp_snarl_records.size() == snarl_index.second+1);
#endif
            temp_index.temp_snarl_records.pop_back();
        } else {
            //This is the child of a chain
            
#ifdef debug_distance_indexing
            cerr << "  Ending new snarl " << temp_index.structure_start_end_as_string(snarl_index)
                 << endl << "    that is a child of " << temp_index.structure_start_end_as_string(temp_snarl_record.parent) << endl;
#endif

#ifdef debug_distance_indexing
            assert(stack.back().first == SnarlDistanceIndex::TEMP_CHAIN);
#endif
            temp_snarl_record.parent = stack.back();
            auto& temp_chain = temp_index.get_chain(stack.back());
            temp_chain.children.emplace_back(snarl_index);
            temp_chain.children.emplace_back(SnarlDistanceIndex::TEMP_NODE, node_id);

        }

        //Record the node itself. This gets done for the start of the chain, and ends of snarls
        SnarlDistanceIndex::temp_record_ref_t node_index = make_pair(SnarlDistanceIndex::TEMP_NODE, node_id);
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index.get_node(node_index);
        temp_node_record.node_id = node_id;
        temp_node_record.node_length = graph->get_length(snarl_end_handle);
        temp_node_record.reversed_in_parent = graph->get_is_reverse(snarl_end_handle);
        temp_node_record.parent = stack.back();
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
        temp_index.components.emplace_back(SnarlDistanceIndex::TEMP_SNARL, temp_index.temp_snarl_records.size());
        temp_index.temp_snarl_records.emplace_back();
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.temp_snarl_records.back();
        temp_snarl_record.is_root_snarl = true;
        temp_snarl_record.parent = make_pair(SnarlDistanceIndex::TEMP_ROOT, 0);


        for (size_t chain_i : root_snarl_indexes) {
            //For each chain component of this root-level snarl
            if (temp_index.root_snarl_components[chain_i].first == SnarlDistanceIndex::TEMP_CHAIN){
                SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = temp_index.get_chain(temp_index.root_snarl_components[chain_i]);
                temp_chain_record.parent = make_pair(SnarlDistanceIndex::TEMP_SNARL, temp_index.temp_snarl_records.size() - 1);
                temp_chain_record.rank_in_parent = temp_snarl_record.children.size();
                temp_chain_record.reversed_in_parent = false;

                temp_snarl_record.children.emplace_back(temp_index.root_snarl_components[chain_i]);
            } else {
#ifdef debug_distance_indexing
                assert(temp_index.root_snarl_components[chain_i].first == SnarlDistanceIndex::TEMP_NODE);
#endif
                SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index.get_node(temp_index.root_snarl_components[chain_i]);
                temp_node_record.parent = make_pair(SnarlDistanceIndex::TEMP_SNARL, temp_index.temp_snarl_records.size() - 1);
                temp_node_record.rank_in_parent = temp_snarl_record.children.size();
                temp_node_record.reversed_in_parent = false;

                temp_snarl_record.children.emplace_back(temp_index.root_snarl_components[chain_i]);
            }
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
    for (int i = temp_index.temp_chain_records.size()-1 ; i >= 0 ; i--) {
        SnarlDistanceIndex::temp_record_ref_t chain_index = make_pair(SnarlDistanceIndex::TEMP_CHAIN, i);
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = temp_index.get_chain(chain_index);
#ifdef debug_distance_indexing
        assert(!temp_chain_record.is_trivial);
        cerr << "  At"  << (temp_chain_record.is_trivial ? " trivial " : "") << "chain " << temp_index.structure_start_end_as_string(chain_index) << endl;
#endif

        //Add the first values for the prefix sum and backwards loop vectors
        temp_chain_record.prefix_sum.emplace_back(0);
        temp_chain_record.max_prefix_sum.emplace_back(0);
        temp_chain_record.backward_loops.emplace_back(std::numeric_limits<size_t>::max());
        temp_chain_record.chain_components.emplace_back(0);


        /*First, go through each of the snarls in the chain in the forward direction and
         * fill in the distances in the snarl. Also fill in the prefix sum and backwards
         * loop vectors here
         */
        size_t curr_component = 0; //which component of the chain are we in
        size_t last_node_length = 0;
        for (size_t chain_child_i = 0 ; chain_child_i < temp_chain_record.children.size() ; chain_child_i++ ){
            const SnarlDistanceIndex::temp_record_ref_t& chain_child_index = temp_chain_record.children[chain_child_i];
            //Go through each of the children in the chain, skipping nodes
            //The snarl may be trivial, in which case don't fill in the distances
#ifdef debug_distance_indexing
            cerr << "    Looking at child " << temp_index.structure_start_end_as_string(chain_child_index) 
                 << " current max prefix sum " << temp_chain_record.max_prefix_sum.back() << endl;
#endif

            if (chain_child_index.first == SnarlDistanceIndex::TEMP_SNARL){
                //This is where all the work gets done. Need to go through the snarl and add
                //all distances, then add distances to the chain that this is in
                //The parent chain will be the last thing in the stack
                SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = 
                        temp_index.get_snarl(chain_child_index);

                //Fill in this snarl's distances
                populate_snarl_index(temp_index, chain_child_index, size_limit, only_top_level_chain_distances, graph);

                bool new_component = temp_snarl_record.min_length == std::numeric_limits<size_t>::max();
                if (new_component){
                    curr_component++;
                }

                //And get the distance values for the end node of the snarl in the chain
                if (new_component) {
                    //If this snarl wasn't start-end connected, then we start 
                    //tracking the distance vectors here

                    //Update the maximum distance
                    temp_index.max_distance = std::max(temp_index.max_distance, temp_chain_record.max_prefix_sum.back());

                    temp_chain_record.prefix_sum.emplace_back(0);
                    temp_chain_record.max_prefix_sum.emplace_back(0);
                    temp_chain_record.backward_loops.emplace_back(temp_snarl_record.distance_end_end);
                    //If the chain is disconnected, the max length is infinite
                    temp_chain_record.max_length =  std::numeric_limits<size_t>::max();
                } else {
                    temp_chain_record.prefix_sum.emplace_back(SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                              temp_chain_record.prefix_sum.back(),
                                                              temp_snarl_record.min_length), 
                                                              temp_snarl_record.start_node_length));
                    temp_chain_record.max_prefix_sum.emplace_back(SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                                   temp_chain_record.max_prefix_sum.back(),
                                                                   temp_snarl_record.max_length), 
                                                                   temp_snarl_record.start_node_length));
                    temp_chain_record.backward_loops.emplace_back(std::min(temp_snarl_record.distance_end_end,
                        SnarlDistanceIndex::sum(temp_chain_record.backward_loops.back()
                        , 2 * (temp_snarl_record.start_node_length + temp_snarl_record.min_length))));
                    temp_chain_record.max_length = SnarlDistanceIndex::sum(temp_chain_record.max_length,
                                                                           temp_snarl_record.max_length);
                }
                temp_chain_record.chain_components.emplace_back(curr_component);
                if (chain_child_i == temp_chain_record.children.size() - 2 && temp_snarl_record.min_length == std::numeric_limits<size_t>::max()) {
                    temp_chain_record.loopable = false;
                }
                last_node_length = 0;
            } else {
                if (last_node_length != 0) {
                    //If this is a node and the last thing was also a node,
                    //then there was a trivial snarl 
                    SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                            temp_index.get_node(chain_child_index);

                    //Check if there is a loop in this node
                    //Snarls get counted as trivial if they contain no nodes but they might still have edges
                    size_t backward_loop = std::numeric_limits<size_t>::max();

                    graph->follow_edges(graph->get_handle(temp_node_record.node_id, !temp_node_record.reversed_in_parent), false, [&](const handle_t& next_handle) {
                        if (graph->get_id(next_handle) == temp_node_record.node_id) {
                            //If there is a loop going backwards (relative to the chain) back to the same node
                            backward_loop = 0;
                        }
                    });

                    temp_chain_record.prefix_sum.emplace_back(SnarlDistanceIndex::sum(temp_chain_record.prefix_sum.back(), last_node_length));
                    temp_chain_record.max_prefix_sum.emplace_back(SnarlDistanceIndex::sum(temp_chain_record.max_prefix_sum.back(), last_node_length));
                    temp_chain_record.backward_loops.emplace_back(std::min(backward_loop,
                        SnarlDistanceIndex::sum(temp_chain_record.backward_loops.back(), 2 * last_node_length)));

                    if (chain_child_i == temp_chain_record.children.size()-1) {
                        //If this is the last node
                        temp_chain_record.loopable=false;
                    }
                    temp_chain_record.chain_components.emplace_back(curr_component);
                }
                last_node_length = temp_index.get_node(chain_child_index).node_length;
                //And update the chains max length
                temp_chain_record.max_length = SnarlDistanceIndex::sum(temp_chain_record.max_length,
                                                                       last_node_length);
            }
        } //Finished walking through chain
        if (temp_chain_record.start_node_id == temp_chain_record.end_node_id && temp_chain_record.chain_components.back() != 0) {
            //If this is a looping, multicomponent chain, the start/end node could end up in separate chain components
            //despite being the same node.
            //Since the first component will always be 0, set the first node's component to be whatever the last
            //component was
            temp_chain_record.chain_components[0] = temp_chain_record.chain_components.back();

        }

        //For a multicomponent chain, the actual minimum length will always be infinite, but since we sometimes need
        //the length of the last component, save that here
        temp_chain_record.min_length = !temp_chain_record.is_trivial && temp_chain_record.start_node_id == temp_chain_record.end_node_id
                        ? temp_chain_record.prefix_sum.back()
                        : SnarlDistanceIndex::sum(temp_chain_record.prefix_sum.back() , temp_chain_record.end_node_length);

#ifdef debug_distance_indexing
        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.backward_loops.size());
        assert(temp_chain_record.prefix_sum.size() == temp_chain_record.chain_components.size());
#endif


        /*Now that we've gone through all the snarls in the chain, fill in the forward loop vector
         * by going through the chain in the backwards direction
         */
        temp_chain_record.forward_loops.resize(temp_chain_record.prefix_sum.size(),
                                               std::numeric_limits<size_t>::max());
        if (temp_chain_record.start_node_id == temp_chain_record.end_node_id && temp_chain_record.children.size() > 1) {

            //If this is a looping chain, then check the first snarl for a loop
            if (temp_chain_record.children.at(1).first == SnarlDistanceIndex::TEMP_SNARL) {
                SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(temp_chain_record.children.at(1));
                temp_chain_record.forward_loops[temp_chain_record.forward_loops.size()-1] = temp_snarl_record.distance_start_start;
            } 
        }

        size_t node_i = temp_chain_record.prefix_sum.size() - 2;
        // We start at the next to last node because we need to look at this record and the next one.
        last_node_length = 0;
        for (int j = (int)temp_chain_record.children.size() - 1 ; j >= 0 ; j--) {
            auto& child = temp_chain_record.children.at(j);
            if (child.first == SnarlDistanceIndex::TEMP_SNARL){
                SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(child);
                if (temp_chain_record.chain_components.at(node_i) != temp_chain_record.chain_components.at(node_i+1) &&
                    temp_chain_record.chain_components.at(node_i+1) != 0){
                    //If this is a new chain component, then add the loop distance from the snarl
                    //If the component of the next node is 0, then we're still in the same component since we're going backwards
                    temp_chain_record.forward_loops.at(node_i) = temp_snarl_record.distance_start_start;
                } else {
                    temp_chain_record.forward_loops.at(node_i) =
                        std::min(SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                    temp_chain_record.forward_loops.at(node_i+1), 
                                    2* temp_snarl_record.min_length),
                                    2*temp_snarl_record.end_node_length), 
                                temp_snarl_record.distance_start_start);
                }
                node_i --;
                last_node_length = 0;
            } else {
                if (last_node_length != 0) {
                    SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = 
                            temp_index.get_node(child);


                    //Check if there is a loop in this node
                    //Snarls get counted as trivial if they contain no nodes but they might still have edges
                    size_t forward_loop = std::numeric_limits<size_t>::max();
                    graph->follow_edges(graph->get_handle(temp_node_record.node_id, temp_node_record.reversed_in_parent), false, [&](const handle_t& next_handle) {
                        if (graph->get_id(next_handle) == temp_node_record.node_id) {
                            //If there is a loop going forward (relative to the chain) back to the same node
                            forward_loop = 0;
                        }
                    });
                    temp_chain_record.forward_loops.at(node_i) = std::min( forward_loop,
                        SnarlDistanceIndex::sum(temp_chain_record.forward_loops.at(node_i+1) , 
                                                 2*last_node_length));
                    node_i--;
                }
                last_node_length = temp_index.get_node(child).node_length;
            }
        }


        //If this is a looping chain, check if the loop distances can be improved by going around the chain

        if (temp_chain_record.start_node_id == temp_chain_record.end_node_id && temp_chain_record.children.size() > 1) {


            //Also check if the reverse loop values would be improved if we went around again

            if (temp_chain_record.backward_loops.back() < temp_chain_record.backward_loops.front()) {
                temp_chain_record.backward_loops[0] = temp_chain_record.backward_loops.back();
                size_t node_i = 1;
                size_t last_node_length = 0;
                for (size_t i = 1 ; i < temp_chain_record.children.size()-1 ; i++ ) {
                    auto& child = temp_chain_record.children.at(i);
                    if (child.first == SnarlDistanceIndex::TEMP_SNARL) {
                        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(child);
                        size_t new_loop_distance = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                      temp_chain_record.backward_loops.at(node_i-1), 
                                                      2*temp_snarl_record.min_length), 
                                                      2*temp_snarl_record.start_node_length); 
                        if (temp_chain_record.chain_components.at(node_i)!= 0 || new_loop_distance >= temp_chain_record.backward_loops.at(node_i)) {
                            //If this is a new chain component or it doesn't improve, stop
                            break;
                        } else {
                            //otherwise record the better distance
                            temp_chain_record.backward_loops.at(node_i) = new_loop_distance;

                        }
                        node_i++;
                        last_node_length = 0;
                    } else {
                        if (last_node_length != 0) {
                            size_t new_loop_distance = SnarlDistanceIndex::sum(temp_chain_record.backward_loops.at(node_i-1), 
                                    2*last_node_length); 
                            size_t old_loop_distance = temp_chain_record.backward_loops.at(node_i);
                            temp_chain_record.backward_loops.at(node_i) = std::min(old_loop_distance,new_loop_distance);
                            node_i++;
                        }
                        last_node_length = temp_index.get_node(child).node_length;
                    }
                }
            }
            if (temp_chain_record.forward_loops.front() < temp_chain_record.forward_loops.back()) {
                //If this is a looping chain and looping improves the forward loops, 
                //then we have to keep going around to update distance

                temp_chain_record.forward_loops.back() = temp_chain_record.forward_loops.front();
                size_t last_node_length = 0;
                node_i = temp_chain_record.prefix_sum.size() - 2;
                for (int j = (int)temp_chain_record.children.size() - 1 ; j >= 0 ; j--) {
                    auto& child = temp_chain_record.children.at(j);
                    if (child.first == SnarlDistanceIndex::TEMP_SNARL){
                        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(child);
                        size_t new_distance = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                temp_chain_record.forward_loops.at(node_i+1), 
                                                2* temp_snarl_record.min_length),
                                                2*temp_snarl_record.end_node_length);
                        if (temp_chain_record.chain_components.at(node_i) != temp_chain_record.chain_components.at(node_i+1) ||
                            new_distance >= temp_chain_record.forward_loops.at(node_i)){
                            //If this is a new component or the distance doesn't improve, stop looking
                            break;
                        } else {
                            //otherwise, update the distance
                            temp_chain_record.forward_loops.at(node_i) = new_distance;
                        }
                        node_i --;
                        last_node_length =0;
                    } else {
                        if (last_node_length != 0) {
                            size_t new_distance = SnarlDistanceIndex::sum(temp_chain_record.forward_loops.at(node_i+1) , 2* last_node_length);
                            size_t old_distance = temp_chain_record.forward_loops.at(node_i);
                            temp_chain_record.forward_loops.at(node_i) = std::min(old_distance, new_distance);
                            node_i--;
                        }
                        last_node_length = temp_index.get_node(child).node_length;
                    }
                } 
            }
        }

        temp_index.max_distance = std::max(temp_index.max_distance, temp_chain_record.max_prefix_sum.back());
        temp_index.max_distance = temp_chain_record.forward_loops.back() == std::numeric_limits<size_t>::max() ? temp_index.max_distance : std::max(temp_index.max_distance, temp_chain_record.forward_loops.back());
        temp_index.max_distance = temp_chain_record.backward_loops.front() == std::numeric_limits<size_t>::max() ? temp_index.max_distance : std::max(temp_index.max_distance, temp_chain_record.backward_loops.front());
        assert(temp_index.max_distance <= 2742664019);

    }

#ifdef debug_distance_indexing
    cerr << "Filling in the distances in root snarls and distances along chains" << endl;
#endif
    for (SnarlDistanceIndex::temp_record_ref_t& component_index : temp_index.components) {
        if (component_index.first == SnarlDistanceIndex::TEMP_SNARL) {
            SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index.get_snarl(component_index);
            populate_snarl_index(temp_index, component_index, size_limit, only_top_level_chain_distances, graph);
            temp_snarl_record.min_length = std::numeric_limits<size_t>::max();
        }
    }
    temp_index.root_structure_count = temp_index.components.size();
#ifdef debug_distance_indexing
    assert(temp_index.components.size() == temp_index.root_structure_count);
    cerr << "Finished temp index with " << temp_index.root_structure_count << " connected components" << endl;
#endif
    return temp_index;
}

/**
 * Populate a row of the distance matrix.
 * Also responsible for filling in min_length, distance_start_start, and distance_start_end on the TemporarySnarlRecord when a distance matrix is used.
 */
static void populate_distance_matrix_row(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const SnarlDistanceIndex::temp_record_ref_t& start_index, const HandleGraph* graph, size_t start_rank, bool is_internal_node, size_t size_limit); 

/** 
 * Fills in required distance matrix rows for each child
 * - Normal snarl: all rows
 * - Oversized snarl: boundaries and tips
 * - size_limit == 0: no distances in index, so no rows
 * - Top-level chain distances only: ??? 
 */
static void populate_distance_matrix_if_needed(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph, size_t size_limit, bool only_top_level_chain_distances); 

/**
 * Does three things:
 * - Builds temp graph that hub labels will be built on
 * - Builds the hub labels
 * - Stores labels in temp_snarl_record
 */
static void populate_hub_labeling(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph);

/**
 * Determine if a snarl is regular or not.
 *
 * A regular snarl is a snarl that consists of only nodes or
 * chains connected to the start and end, without any connections between
 * multiple children, or any way to turn around. There may be an edge directly
 * across.
 *
 * A simple snarl is always regular.
 */
bool check_regularity(const SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, const SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph);

// ---------------------------------------------------------------------------
// Phase helpers for populate_snarl_index (all file-static)
// ---------------------------------------------------------------------------

// Step 1: Walk up the snarl tree from curr_index until we find the direct
// child of ancestor_snarl_index that contains curr_index.
static SnarlDistanceIndex::temp_record_ref_t ancestor_of_node_in_snarl(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        SnarlDistanceIndex::temp_record_ref_t curr_index,
        SnarlDistanceIndex::temp_record_ref_t ancestor_snarl_index) {

    const auto& snarl = temp_index.get_snarl(ancestor_snarl_index);
    if (curr_index.second == snarl.start_node_id ||
        curr_index.second == snarl.end_node_id) {
        return curr_index;
    }

    SnarlDistanceIndex::temp_record_ref_t parent_index = temp_index.get_node(curr_index).parent;
    while (parent_index != ancestor_snarl_index) {
        curr_index = parent_index;
        parent_index = parent_index.first == SnarlDistanceIndex::TEMP_SNARL
                        ? temp_index.get_snarl(parent_index).parent
                        : temp_index.get_chain(parent_index).parent;
#ifdef debug_distance_indexing
        assert(parent_index.first != SnarlDistanceIndex::TEMP_ROOT);
#endif
    }
    return curr_index;
}

// Step 2a: Return the handle pointing out from child_index in the given
// traversal direction (reversed=false → forward/end side; reversed=true →
// backward/start side).
static handle_t child_boundary_handle(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::temp_record_ref_t& child_index,
        bool reversed,
        const HandleGraph* graph) {

    if (child_index.first == SnarlDistanceIndex::TEMP_NODE) {
        return graph->get_handle(child_index.second, reversed);
    } else if (reversed) {
        return graph->get_handle(temp_index.get_chain(child_index).start_node_id,
                                 !temp_index.get_chain(child_index).start_node_rev);
    } else {
        return graph->get_handle(temp_index.get_chain(child_index).end_node_id,
                                  temp_index.get_chain(child_index).end_node_rev);
    }
}

// Step 2b: Determine the traversal direction of child_index when entered via
// graph_handle.
static bool child_side_reversed(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::temp_record_ref_t& child_index,
        handle_t graph_handle,
        const HandleGraph* graph) {

    if (child_index.first == SnarlDistanceIndex::TEMP_NODE ||
        temp_index.get_chain(child_index).is_trivial) {
        return graph->get_is_reverse(graph_handle);
    }
    return graph->get_id(graph_handle) == temp_index.get_chain(child_index).end_node_id;
}

// ---------------------------------------------------------------------------
// SnarlChildGraph implementation
// ---------------------------------------------------------------------------

SnarlChildGraph::SnarlChildGraph(
        TempIndex& temp_index,
        temp_record_ref_t snarl_index,
        std::span<const temp_record_ref_t> children,
        const handlegraph::HandleGraph* graph)
    : temp_index_(temp_index)
    , snarl_index_(snarl_index)
    , children_(children)
    , graph_(graph) {}

std::span<const SnarlChildGraph::temp_record_ref_t>
SnarlChildGraph::children() const noexcept {
    return children_;
}

std::pair<SnarlChildGraph::temp_record_ref_t, bool>
SnarlChildGraph::boundary(bool start) const {
    const auto& snarl = temp_index_.get_snarl(snarl_index_);
    if (start) {
        return {{SnarlDistanceIndex::TEMP_NODE, snarl.start_node_id}, snarl.start_node_rev};
    }
    return {{SnarlDistanceIndex::TEMP_NODE, snarl.end_node_id}, snarl.end_node_rev};
}

void SnarlChildGraph::for_each_outgoing(
        temp_record_ref_t child,
        bool go_left,
        const std::function<void(temp_record_ref_t neighbor,
                                 bool neighbor_rev,
                                 size_t edge_distance,
                                 handlegraph::nid_t arriving_node_id)>& callback) const {

    handle_t out_handle = child_boundary_handle(temp_index_, child, go_left, graph_);

    graph_->follow_edges(out_handle, false, [&](const handle_t& next_handle) {
        handlegraph::nid_t arriving_nid = graph_->get_id(next_handle);
        temp_record_ref_t next_node = {SnarlDistanceIndex::TEMP_NODE, arriving_nid};
        temp_record_ref_t neighbor  = ancestor_of_node_in_snarl(temp_index_, next_node, snarl_index_);
        bool neighbor_rev           = child_side_reversed(temp_index_, neighbor, next_handle, graph_);

        size_t edge_distance;
        if (neighbor.first == SnarlDistanceIndex::TEMP_NODE) {
            edge_distance = graph_->get_length(next_handle);
        } else {
            const auto& chain = temp_index_.get_chain(neighbor);
            edge_distance = chain.min_length;
            if (chain.chain_components.back() != 0) {
                edge_distance = std::numeric_limits<size_t>::max();
            }
        }

        callback(neighbor, neighbor_rev, edge_distance, arriving_nid);
        return true;
    });
}

// Phase 1: Mark tip nodes and set is_simple=false if any tip is found.
static void identify_tips(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
        const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children,
        const HandleGraph* graph) {

    for (const auto& child : all_children) {
        if (child.first != SnarlDistanceIndex::TEMP_NODE
            || (child.second != temp_snarl_record.start_node_id
                && child.second != temp_snarl_record.end_node_id)) {
            bool is_node = (child.first == SnarlDistanceIndex::TEMP_NODE);
            nid_t node_id = is_node ? child.second
                                    : temp_index.temp_chain_records.at(child.second).end_node_id;
            size_t rank = is_node ? temp_index.temp_node_records.at(child.second - temp_index.min_node_id).rank_in_parent
                                  : temp_index.temp_chain_records.at(child.second).rank_in_parent;
            bool is_reverse = is_node ? false
                                      : temp_index.temp_chain_records.at(child.second).end_node_rev;
            rank -= 2;

            bool has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, is_reverse), false, [&](const handle_t next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_index.temp_node_records.at(node_id - temp_index.min_node_id).is_tip = true;
                temp_snarl_record.tippy_child_ranks.emplace(rank, false);
                temp_snarl_record.is_simple = false;
            }
            node_id = is_node ? child.second
                              : temp_index.temp_chain_records.at(child.second).start_node_id;
            is_reverse = is_node ? true
                                 : !temp_index.temp_chain_records.at(child.second).start_node_rev;
            has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, is_reverse), false, [&](const handle_t next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_index.temp_node_records.at(node_id - temp_index.min_node_id).is_tip = true;
                temp_snarl_record.tippy_child_ranks.emplace(rank, true);
                temp_snarl_record.is_simple = false;
            }
        }
    }
}

// Phase 2a: BFS topological sort of children; returns new-to-old rank mapping.
// TODO: For non-DAGs this sort will end up arbitrary. That doesn't matter
//       since the only consumer of ranks (ziptrees) expects arbitrary ranks.
static vector<size_t> topo_sort_children(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::temp_record_ref_t& snarl_index,
        const SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
        std::span<const SnarlDistanceIndex::temp_record_ref_t> all_children,
        const HandleGraph* graph) {

    SnarlChildGraph child_graph(temp_index, snarl_index, all_children, graph);

    vector<size_t> topological_sort_order;
    topological_sort_order.reserve(all_children.size());

    unordered_set<size_t> visited_ranks;
    visited_ranks.reserve(all_children.size());

    vector<pair<size_t, bool>> source_nodes;

    // Add tips as sources. Tips push before the start sentinel so sentinel pops first (LIFO).
    for (const auto& tip : temp_snarl_record.tippy_child_ranks) {
        source_nodes.emplace_back(tip.first, !tip.second);
    }
    // Start node dummy rank is max() — pops first as LIFO sentinel.
    source_nodes.emplace_back(std::numeric_limits<size_t>::max(), false);

    while (!source_nodes.empty()) {
        pair<size_t, bool> current_child_index = source_nodes.back();
        source_nodes.pop_back();

        if (visited_ranks.count(current_child_index.first) != 0) {
            // Revisiting a source means we hit a loop; abort with arbitrary ranks.
            break;
        }
        if (current_child_index.first != std::numeric_limits<size_t>::max()) {
            topological_sort_order.emplace_back(current_child_index.first);
        }
        visited_ranks.emplace(current_child_index.first);

        // Determine which child (or sentinel start boundary) to follow edges from.
        // For the sentinel, use the snarl's start boundary node in its stored orientation,
        // which produces the same handle as topological_sort_start in the original code.
        SnarlDistanceIndex::temp_record_ref_t current_ref;
        bool go_left;
        if (current_child_index.first == std::numeric_limits<size_t>::max()) {
            current_ref = {SnarlDistanceIndex::TEMP_NODE, temp_snarl_record.start_node_id};
            go_left     = temp_snarl_record.start_node_rev;
        } else {
            current_ref = all_children[current_child_index.first];
            go_left     = current_child_index.second;
        }

        child_graph.for_each_outgoing(current_ref, go_left, [&](
                SnarlDistanceIndex::temp_record_ref_t neighbor,
                bool neighbor_rev,
                size_t /*edge_distance*/,
                handlegraph::nid_t /*arriving_nid*/) {
#ifdef debug_distance_indexing
            cerr << "Following forward edges to "
                 << temp_index.structure_start_end_as_string(neighbor) << endl;
#endif
            // Skip snarl boundaries.
            if (neighbor.first == SnarlDistanceIndex::TEMP_NODE &&
                (neighbor.second == temp_snarl_record.start_node_id ||
                 neighbor.second == temp_snarl_record.end_node_id)) {
                return;
            }
            size_t next_rank = neighbor.first == SnarlDistanceIndex::TEMP_NODE
                        ? temp_index.get_node(neighbor).rank_in_parent
                        : temp_index.get_chain(neighbor).rank_in_parent;
            assert(next_rank >= 2);
            next_rank -= 2;
            assert(all_children[next_rank] == neighbor);
            if (visited_ranks.count(next_rank) != 0) {
                return;
            }

            // Check if neighbor is a topological source (no unvisited predecessors).
            bool is_source = true;
            child_graph.for_each_outgoing(neighbor, !neighbor_rev, [&](
                    SnarlDistanceIndex::temp_record_ref_t incoming,
                    bool /*incoming_rev*/,
                    size_t /*edge_distance*/,
                    handlegraph::nid_t /*arriving_nid*/) {
#ifdef debug_distance_indexing
                cerr << "Getting backwards edge from "
                     << temp_index.structure_start_end_as_string(incoming) << endl;
#endif
                if (incoming.first == SnarlDistanceIndex::TEMP_NODE &&
                    (incoming.second == temp_snarl_record.start_node_id ||
                     incoming.second == temp_snarl_record.end_node_id)) {
                    return;
                }
                size_t incoming_rank = incoming.first == SnarlDistanceIndex::TEMP_NODE
                            ? temp_index.get_node(incoming).rank_in_parent
                            : temp_index.get_chain(incoming).rank_in_parent;
                assert(incoming_rank >= 2);
                incoming_rank -= 2;
                if (visited_ranks.count(incoming_rank) == 0) {
                    is_source = false;
                }
            });
            if (is_source) {
                source_nodes.emplace_back(next_rank, neighbor_rev);
            }
        });
    }

    // Non-DAG fallback: append any ranks not yet visited in arbitrary order.
    vector<bool> check_ranks(all_children.size(), false);
    for (size_t x : topological_sort_order) {
        check_ranks[x] = true;
    }
    for (size_t i = 0; i < check_ranks.size(); i++) {
        if (!check_ranks[i]) {
            topological_sort_order.emplace_back(i);
        }
    }
    assert(topological_sort_order.size() == all_children.size());
    return topological_sort_order;
}

// Phase 2b: Apply the topo-sort permutation: update rank_in_parent for every
// child and rebuild tippy_child_ranks with new ranks.
static void apply_topo_permutation(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
        std::span<const SnarlDistanceIndex::temp_record_ref_t> all_children,
        const vector<size_t>& new_to_old) {

    auto old_tippy_ranks = temp_snarl_record.tippy_child_ranks;
    temp_snarl_record.tippy_child_ranks.clear();
    for (size_t new_rank = 0; new_rank < new_to_old.size(); new_rank++) {
        size_t old_rank = new_to_old[new_rank];
        if (all_children[old_rank].first == SnarlDistanceIndex::TEMP_NODE) {
            temp_index.get_node(all_children[old_rank]).rank_in_parent = new_rank + 2;
        } else {
            temp_index.get_chain(all_children[old_rank]).rank_in_parent = new_rank + 2;
        }
        const auto& old_is_tip = old_tippy_ranks.find(old_rank);
        if (old_is_tip != old_tippy_ranks.end()) {
            temp_snarl_record.tippy_child_ranks.emplace(new_rank, old_is_tip->second);
        }
    }
}

// Phase 3: Compute snarl distances (normal or oversized hub-label path).
// Appends boundary nodes to all_children (unless is_root_snarl).
static void compute_snarl_distances(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::temp_record_ref_t& snarl_index,
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
        vector<SnarlDistanceIndex::temp_record_ref_t>& all_children,
        const HandleGraph* graph,
        size_t size_limit,
        bool only_top_level_chain_distances) {

    if (!temp_snarl_record.is_root_snarl) {
        all_children.emplace_back(SnarlDistanceIndex::TEMP_NODE, temp_snarl_record.start_node_id);
        all_children.emplace_back(SnarlDistanceIndex::TEMP_NODE, temp_snarl_record.end_node_id);
    }

    if (size_limit != 0 && temp_snarl_record.node_count > size_limit) {
        temp_index.most_oversized_snarl_size = std::max(temp_index.most_oversized_snarl_size, temp_snarl_record.node_count);
        temp_index.use_oversized_snarls = true;
        temp_snarl_record.is_simple = false;
        populate_hub_labeling(temp_index, snarl_index, temp_snarl_record, all_children, graph);

        // Query hub labeling for connectivity distances (excluding boundary lengths).
        // Start is always child rank 0 forward, end is always child rank 1 forward.
        temp_snarl_record.min_length = promote_distance<size_t>(hhl_query(temp_snarl_record.hub_labels.begin(), bgid(0, false, true), bgid(1, false, false)));
        temp_snarl_record.distance_start_start = promote_distance<size_t>(hhl_query(temp_snarl_record.hub_labels.begin(), bgid(0, false, true), bgid(0, true, false)));
        temp_snarl_record.distance_end_end = promote_distance<size_t>(hhl_query(temp_snarl_record.hub_labels.begin(), bgid(1, true, true), bgid(1, false, false)));
        // TODO: Should this be here or should it be part of populate_hub_labeling()? Or its own function?
    } else {
        if (size_limit == 0 || only_top_level_chain_distances) {
            temp_snarl_record.include_distances = false;
        }
        // Also fills in min_length, distance_start_start, distance_start_end, sets is_simple=false if not simple.
        populate_distance_matrix_if_needed(temp_index, snarl_index, temp_snarl_record, all_children, graph, size_limit, only_top_level_chain_distances);
    }
}

// Phase 4: For simple snarls, record child node orientations.
// IMPORTANT: iterates temp_snarl_record.children[0..node_count), NOT all_children —
// boundary nodes appended in Phase 3 must not be included here.
static void mark_simple_snarl_orientations(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record) {

    for (size_t i = 0; i < temp_snarl_record.node_count; i++) {
        const SnarlDistanceIndex::temp_record_ref_t& child_index = temp_snarl_record.children[i];
#ifdef debug_distance_indexing
        assert(child_index.first == SnarlDistanceIndex::TEMP_NODE);
#endif
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record =
            temp_index.get_node(child_index);
        temp_node_record.reversed_in_parent =
            temp_node_record.distance_left_start == std::numeric_limits<size_t>::max();
    }
}

// Phase 5: Regularity check and index-size accounting.
static void finalize_snarl_record(
        SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
        const SnarlDistanceIndex::temp_record_ref_t& snarl_index,
        SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record,
        const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children,
        const HandleGraph* graph) {

    temp_snarl_record.is_regular = check_regularity(temp_index, snarl_index, temp_snarl_record, all_children, graph);

    temp_index.max_index_size += temp_snarl_record.get_max_record_length();
    if (temp_snarl_record.is_simple) {
        temp_index.max_index_size -= (temp_snarl_record.children.size() *
            SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryNodeRecord::get_max_record_length());
    }
    temp_index.max_bits = std::max(temp_index.max_bits,
        22 + SnarlDistanceIndex::bit_width(temp_snarl_record.children.size()));
}

/**
 * Fill in the snarl index.
 * The index will already know its boundaries and everything knows their relationships in the
 * snarl tree. This needs to fill in the distances and the ranks of children in the snarl
 * The rank of a child is arbitrary, except that the start node will always be 0 and the end node
 * will always be the node count+1 (since node count doesn't count the boundary nodes)
 */
void populate_snarl_index(
                SnarlDistanceIndex::TemporaryDistanceIndex& temp_index,
                SnarlDistanceIndex::temp_record_ref_t snarl_index, size_t size_limit,
                bool only_top_level_chain_distances, const HandleGraph* graph) {
#ifdef debug_distance_indexing
    cerr << "Getting the distances for snarl " << temp_index.structure_start_end_as_string(snarl_index) << endl;
    assert(snarl_index.first == SnarlDistanceIndex::TEMP_SNARL);
#endif
    auto& temp_snarl_record = temp_index.get_snarl(snarl_index);
    temp_snarl_record.is_simple = true;

    vector<SnarlDistanceIndex::temp_record_ref_t> all_children = temp_snarl_record.children;

    identify_tips(temp_index, temp_snarl_record, all_children, graph);

    if (!temp_snarl_record.is_root_snarl) {
        auto new_to_old = topo_sort_children(temp_index, snarl_index, temp_snarl_record, all_children, graph);
        apply_topo_permutation(temp_index, temp_snarl_record, all_children, new_to_old);
    }

    compute_snarl_distances(temp_index, snarl_index, temp_snarl_record, all_children, graph,
                            size_limit, only_top_level_chain_distances);

#ifdef debug_distance_indexing
    cerr << "snarl " << temp_index.structure_start_end_as_string(snarl_index) << " is_simple: " << temp_snarl_record.is_simple << endl;
#endif

    if (temp_snarl_record.is_simple) {
        mark_simple_snarl_orientations(temp_index, temp_snarl_record);
    }

    finalize_snarl_record(temp_index, snarl_index, temp_snarl_record, all_children, graph);
}


void populate_hub_labeling(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph) {
  CHOverlay ov = make_boost_graph(temp_index, snarl_index, temp_snarl_record, all_children, graph);

#ifdef debug_hub_label_build
  // Dump CHOverlay graph to stderr for debugging
  std::cerr << "=== CHOverlay Graph Dump ===" << std::endl;
  std::cerr << ov << std::endl;
  std::cerr << "=== End CHOverlay Dump ===" << std::endl;
#endif

  make_contraction_hierarchy(ov);

  vector<vector<HubRecord>> labels; labels.resize(num_vertices(ov));
  vector<vector<HubRecord>> labels_rev; labels_rev.resize(num_vertices(ov)); 
  create_labels(labels, labels_rev, ov);
#ifdef debug_hub_label_storage
  std::cerr << "Hub labels unpacked:" << std::endl;
  for (const auto& node_list : {labels, labels_rev}) {
    std::cerr << "Labels for all nodes:" << std::endl;
    for (size_t i = 0; i < node_list.size(); i++) {
        std::cerr << "\tLabels for rank " << i << ":" << std::endl;
        for (const HubRecord& label : node_list[i]) {
            std::cerr << "\t\tHub: " << label.hub << " Dist: " << label.dist << std::endl; 
        }
    }
  }
#endif
  
  // Put labels in temp_snarl_record
  temp_snarl_record.hub_labels = pack_labels(labels, labels_rev);
#ifdef debug_hub_label_storage
  std::cerr << "Hub labels as packed: ";
  for (size_t i = 0; i < temp_snarl_record.hub_labels.size(); i++) {
    if (i > 0) {
        std::cerr << " | ";
    }
    std::cerr << temp_snarl_record.hub_labels[i];
  }
  std::cerr << std::endl;
#endif
}

void populate_distance_matrix_if_needed(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph, size_t size_limit, bool only_top_level_chain_distances) {
    if (size_limit != 0 && !only_top_level_chain_distances) { 
      //If we are saving distances
      //Reserve enough space to store all possible distances
      temp_snarl_record.distances.reserve( temp_snarl_record.node_count > size_limit
              ? temp_snarl_record.node_count * 2
              : temp_snarl_record.node_count * temp_snarl_record.node_count);
    } else {
      temp_snarl_record.include_distances = false;
    }
    for (auto it = all_children.rbegin(); it != all_children.rend(); ++it) {
        // Visit all the children in reverse order
        const SnarlDistanceIndex::temp_record_ref_t& start_index = *it;

        bool is_internal_node = false;

        if ((start_index.first == SnarlDistanceIndex::TEMP_NODE 
             && start_index.second != temp_snarl_record.start_node_id 
             && start_index.second != temp_snarl_record.end_node_id) 
            || 
            (start_index.first == SnarlDistanceIndex::TEMP_CHAIN && temp_index.get_chain(start_index).is_trivial)) {
            // If this is an internal node
            is_internal_node = true;
            nid_t node_id = start_index.first == SnarlDistanceIndex::TEMP_NODE ? start_index.second : temp_index.get_chain(start_index).start_node_id;
            SnarlDistanceIndex::temp_record_ref_t node_index {SnarlDistanceIndex::TEMP_NODE, node_id};
            size_t rank = start_index.first == SnarlDistanceIndex::TEMP_NODE ? temp_index.get_node(start_index).rank_in_parent
                                                          : temp_index.get_chain(start_index).rank_in_parent;

            bool has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, false), false, [&](const handle_t& next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_index.get_node(node_index).is_tip = true;
                temp_snarl_record.tippy_child_ranks.emplace(rank, false);
                temp_snarl_record.is_simple=false; //It is a tip so this isn't simple snarl
            }
            has_edges = false;
            graph->follow_edges(graph->get_handle(node_id, true), false, [&](const handle_t& next_handle) {
                has_edges = true;
            });
            if (!has_edges) {
                temp_index.get_node(node_index).is_tip = true;
                temp_snarl_record.tippy_child_ranks.emplace(rank, true);
                temp_snarl_record.is_simple=false; //It is a tip so this isn't simple snarl
            }
        } else if (start_index.first == SnarlDistanceIndex::TEMP_CHAIN && !temp_index.get_chain(start_index).is_trivial) {
            // If this is an internal chain, then it isn't a simple snarl
            temp_snarl_record.is_simple=false;
        }

        bool start_is_tip = start_index.first == SnarlDistanceIndex::TEMP_NODE 
                      ? temp_index.get_node(start_index).is_tip 
                      : temp_index.get_chain(start_index).is_tip;

        size_t start_rank = start_index.first == SnarlDistanceIndex::TEMP_NODE 
                ? temp_index.get_node(start_index).rank_in_parent
                : temp_index.get_chain(start_index).rank_in_parent;


        if (start_index.first == SnarlDistanceIndex::TEMP_NODE && start_index.second == temp_snarl_record.start_node_id) {
            start_rank = 0;
        } else if (start_index.first == SnarlDistanceIndex::TEMP_NODE && start_index.second == temp_snarl_record.end_node_id) {
            start_rank = 1;
        } //TODO:
          //else {
          //  assert(start_rank != 0 && start_rank != 1);
          //} 

        //traversal start is not a tip or a boundary node
        bool start_normal_child = (!start_is_tip && start_rank != 0 && start_rank != 1);
 
        if ( (temp_snarl_record.node_count > size_limit || size_limit == 0 || only_top_level_chain_distances) && (temp_snarl_record.is_root_snarl || start_normal_child)) {
            //If we don't care about internal distances, and we also are not at a boundary or tip
            //TODO: Why do we care about tips specifically?
            continue;
        }
        //getting here means snarl is not oversized
        //fill in all distances for a row
        populate_distance_matrix_row(temp_index, snarl_index, temp_snarl_record, start_index, graph, start_rank, is_internal_node, size_limit);   
    }                                                                                                                    
}      
      
    
                        
void populate_distance_matrix_row(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const SnarlDistanceIndex::temp_record_ref_t& start_index, const HandleGraph* graph, size_t start_rank, bool is_internal_node, size_t size_limit) {

    // SnarlChildGraph encapsulates the follow_edges + ancestor-resolution kernel.
    // Pass an empty children span — for_each_outgoing doesn't use it.
    SnarlChildGraph child_graph(temp_index, snarl_index,
                                std::span<const SnarlDistanceIndex::temp_record_ref_t>{}, graph);

    //Start from either direction for all nodes, but only going in for start and end
    vector<bool> directions;
    if (start_index.first == SnarlDistanceIndex::TEMP_NODE && start_index.second == temp_snarl_record.start_node_id) {
        directions.emplace_back(temp_snarl_record.start_node_rev);
    } else if (start_index.first == SnarlDistanceIndex::TEMP_NODE && start_index.second == temp_snarl_record.end_node_id){
        directions.emplace_back(!temp_snarl_record.end_node_rev);
    } else {
        directions.emplace_back(true);
        directions.emplace_back(false);
    }
    for (bool start_rev : directions) {
        //Start a dijkstra traversal from start_index going in the direction indicated by start_rev
        //Record the distances to each node (child of the snarl) found
        size_t reachable_node_count = 0; //How many nodes can we reach from this node side?

#ifdef debug_distance_indexing
        cerr << "  Starting from child " << temp_index.structure_start_end_as_string(start_index)
             << " going " << (start_rev ? "rev" : "fd") << endl;
#endif

        //Define a NetgraphNode as the value for the priority queue:
        // <distance, <<type of node, index into temp_node/chain_records>, direction>
        using NetgraphNode = pair<size_t, pair<SnarlDistanceIndex::temp_record_ref_t, bool>>;
        auto cmp = [] (const NetgraphNode a, const NetgraphNode b) {
            return a.first > b.first;
        };

        //The priority queue of the next nodes to visit, ordered by the distance
        std::priority_queue<NetgraphNode, vector<NetgraphNode>, decltype(cmp)> queue(cmp);
        //The nodes we've already visited
        unordered_set<pair<SnarlDistanceIndex::temp_record_ref_t, bool>> visited_nodes;
        visited_nodes.reserve(temp_snarl_record.node_count * 2);

        //Start from the current start node
        queue.push(make_pair(0, make_pair(start_index, start_rev)));

        while (!queue.empty()) {

            //Get the current node from the queue and pop it out of the queue
            size_t current_distance = queue.top().first;
            SnarlDistanceIndex::temp_record_ref_t current_index = queue.top().second.first;
            bool current_rev = queue.top().second.second;
            if (visited_nodes.count(queue.top().second)) {
                queue.pop();
                continue;
            }
            visited_nodes.emplace(queue.top().second);
            queue.pop();

            // Pre-compute the current child's outgoing node ID and "other side" node ID
            // for is_simple detection (mirrors the original current_end_handle checks).
            nid_t current_end_nid, current_other_side_nid;
            if (current_index.first == SnarlDistanceIndex::TEMP_NODE) {
                current_end_nid        = current_index.second;
                current_other_side_nid = current_index.second;
            } else {
                const auto& ccr = temp_index.get_chain(current_index);
                current_end_nid        = current_rev ? ccr.start_node_id : ccr.end_node_id;
                current_other_side_nid = current_rev ? ccr.end_node_id   : ccr.start_node_id;
            }

#ifdef debug_distance_indexing
            cerr << "    at child " << temp_index.structure_start_end_as_string(current_index) << " going "
                 << (current_rev ? "rev" : "fd") << " outgoing from node " << current_end_nid << endl;
#endif
            child_graph.for_each_outgoing(current_index, current_rev, [&](
                    SnarlDistanceIndex::temp_record_ref_t next_index,
                    bool next_rev,
                    size_t edge_distance,
                    handlegraph::nid_t arriving_nid) {
#ifdef debug_distance_indexing
                cerr << "      see edge " << current_end_nid
                     << " -> " << arriving_nid << endl;
#endif

                // Detect loops that invalidate simplicity (preserves original checks verbatim).
                if (arriving_nid == current_end_nid || arriving_nid == current_other_side_nid) {
                    //If this loops onto the same node side, or to the other end of the same child
                    temp_snarl_record.is_simple = false;
                } else if (!temp_snarl_record.is_root_snarl && start_rank == 0 &&
                           current_index != start_index &&
                           !(next_index.first == SnarlDistanceIndex::TEMP_NODE &&
                             next_index.second == temp_snarl_record.end_node_id)) {
                    //If the starting point of this traversal was the start of the snarl, the current starting point is not the start node,
                    //and we found another child, then this is not a simple snarl
                    temp_snarl_record.is_simple = false;
                } else if (!temp_snarl_record.is_root_snarl && start_rank == 1 &&
                           current_index != start_index &&
                           !(next_index.first == SnarlDistanceIndex::TEMP_NODE &&
                             next_index.second == temp_snarl_record.start_node_id)) {
                    //If the starting point of this traversal was the end of the snarl, the current starting point is not the end node,
                    //and we found another child, then this is not a simple snarl
                    temp_snarl_record.is_simple = false;
                }

                reachable_node_count++;

                bool next_is_tip = start_index.first == SnarlDistanceIndex::TEMP_NODE
                          ? temp_index.get_node(start_index).is_tip
                          : temp_index.get_chain(start_index).is_tip;

                //The rank and orientation of next in the snarl
                size_t next_rank;
                if (next_index.first == SnarlDistanceIndex::TEMP_NODE && next_index.second == temp_snarl_record.start_node_id) {
#ifdef debug_distance_indexing
                    std::cerr << "        edge arrived at start" << std::endl;
#endif
                    next_rank = 0;
                } else if (next_index.first == SnarlDistanceIndex::TEMP_NODE && next_index.second == temp_snarl_record.end_node_id) {
#ifdef debug_distance_indexing
                    std::cerr << "        edge arrived at end" << std::endl;
#endif
                    next_rank = 1;
                } else {
                    next_rank = next_index.first == SnarlDistanceIndex::TEMP_NODE
                            ? temp_index.get_node(next_index).rank_in_parent
                            : temp_index.get_chain(next_index).rank_in_parent;
                    //If the next thing wasn't a boundary node and this was an internal node, then it isn't a simple snarl
                    if (is_internal_node) {
                        temp_snarl_record.is_simple = false;
                    }
                }//TODO: This won't be true of root snarls
                  //else {
                  //  assert(next_rank != 0 && next_rank != 1);
                  //}

                /**Record the distance **/
                bool start_is_boundary = !temp_snarl_record.is_root_snarl && (start_rank == 0 || start_rank == 1);
                bool next_is_boundary = !temp_snarl_record.is_root_snarl && (next_rank == 0 || next_rank == 1);

                pair<size_t, bool> start = start_is_boundary
                    ? make_pair(start_rank, false) : make_pair(start_rank, !start_rev);
                pair<size_t, bool> next = next_is_boundary
                    ? make_pair(next_rank, false) : make_pair(next_rank, next_rev);

                if (size_limit == 0 && start_is_boundary && next_is_boundary) {
                    // If not measuring distances, we need to use
                    // distance_start_start and distance_end_end as
                    // connectivity flags so we can still detect reversals
                    // within chains and recognize regular snarls.
                    if (start_rank == 0 && next_rank == 0) {
                        temp_snarl_record.distance_start_start = 0;
#ifdef debug_distance_indexing
                        cerr << "        set loop indicator start start distance " << temp_snarl_record.distance_start_start << endl;
#endif
                    } else if (start_rank == 1 && next_rank == 1) {
                        temp_snarl_record.distance_end_end = 0;
#ifdef debug_distance_indexing
                        cerr << "        set loop indicator end end distance " << temp_snarl_record.distance_start_start << endl;
#endif
                    }
                } else if (size_limit != 0 &&
                    (temp_snarl_record.node_count <= size_limit || start_is_boundary || next_is_boundary)) {
                    //If the snarl is too big, then we don't record distances between internal nodes
                    //If we are looking at all distances or we are looking at boundaries
                    bool added_new_distance = false;

                    //Set the distance
                    if (start_is_boundary && next_is_boundary) {
                        //If it is between bounds of the snarl, then the snarl stores it
                        if (start_rank == 0 && next_rank == 0 &&
                            temp_snarl_record.distance_start_start == std::numeric_limits<size_t>::max()) {
                            temp_snarl_record.distance_start_start = current_distance;
#ifdef debug_distance_indexing
                            cerr << "        set start start distance " << temp_snarl_record.distance_start_start << endl;
#endif
                            added_new_distance = true;
                        } else if (start_rank == 1 && next_rank == 1 &&
                                   temp_snarl_record.distance_end_end == std::numeric_limits<size_t>::max()) {
                            temp_snarl_record.distance_end_end = current_distance;
#ifdef debug_distance_indexing
                            cerr << "        set end end distance " << temp_snarl_record.distance_start_start << endl;
#endif
                            added_new_distance = true;
                        } else if (((start_rank == 0 && next_rank == 1) || (start_rank == 1 && next_rank == 0))
                                    && temp_snarl_record.min_length == std::numeric_limits<size_t>::max()){
                            temp_snarl_record.min_length = current_distance;
                            added_new_distance = true;

                        }
                    } else if (start_is_boundary){
                        //If start is a boundary node
                        if (next_index.first == SnarlDistanceIndex::TEMP_NODE) {
                            //Next is a node
                            auto& temp_node_record = temp_index.get_node(next_index);
                            if (start_rank == 0 && !next_rev &&
                                    temp_node_record.distance_left_start == std::numeric_limits<size_t>::max()) {
                                temp_node_record.distance_left_start = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 0 && next_rev &&
                                    temp_node_record.distance_right_start == std::numeric_limits<size_t>::max()) {
                                temp_node_record.distance_right_start = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 1 && !next_rev &&
                                    temp_node_record.distance_left_end == std::numeric_limits<size_t>::max()) {
                                temp_node_record.distance_left_end = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 1 && next_rev &&
                                    temp_node_record.distance_right_end == std::numeric_limits<size_t>::max()) {
                                temp_node_record.distance_right_end = current_distance;
                                added_new_distance = true;
                            }
                        }  else {
                            //Next is a chain
                            auto& temp_chain_record = temp_index.get_chain(next_index);
                            if (start_rank == 0 && !next_rev &&
                                    temp_chain_record.distance_left_start == std::numeric_limits<size_t>::max()) {
                                temp_chain_record.distance_left_start = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 0 && next_rev &&
                                    temp_chain_record.distance_right_start == std::numeric_limits<size_t>::max()) {
                                temp_chain_record.distance_right_start = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 1 && !next_rev &&
                                    temp_chain_record.distance_left_end == std::numeric_limits<size_t>::max()) {
                                temp_chain_record.distance_left_end = current_distance;
                                added_new_distance = true;
                            } else if (start_rank == 1 && next_rev &&
                                    temp_chain_record.distance_right_end == std::numeric_limits<size_t>::max()) {
                                temp_chain_record.distance_right_end = current_distance;
                                added_new_distance = true;
                            }
                        }
                    } else if (!next_is_boundary && !temp_snarl_record.distances.count(make_pair(start, next))) {
                        //Otherwise the snarl stores it in its distance
                        //If the distance isn't from an internal node to a bound and we haven't stored the distance yet

                        temp_snarl_record.distances[make_pair(start, next)] = current_distance;
                        added_new_distance = true;
#ifdef debug_distance_indexing
                        cerr << "           Adding distance between ranks " << start.first << " " << start.second << " and " << next.first << " " << next.second << ": " << current_distance << endl;
#endif
                    }
                    if (added_new_distance) {
                        temp_snarl_record.max_distance = std::max(temp_snarl_record.max_distance, current_distance);
                    }
                }


                /**Add the next node to the priority queue**/

                if (visited_nodes.count(make_pair(next_index, next_rev)) == 0 &&
                    arriving_nid != temp_snarl_record.start_node_id &&
                    arriving_nid != temp_snarl_record.end_node_id
                    ) {
                    //If this isn't leaving the snarl,
                    //then add the next node to the queue, along with the distance to traverse it
                    // edge_distance already encodes chain.min_length (∞ if disconnected).
                    if (edge_distance != std::numeric_limits<size_t>::max()) {
                        queue.push(make_pair(SnarlDistanceIndex::sum(current_distance, edge_distance),
                                       make_pair(next_index, next_rev)));
                    }
                }
                if (next_index.first == SnarlDistanceIndex::TEMP_CHAIN) {
                    size_t loop_distance = next_rev ? temp_index.get_chain(next_index).backward_loops.back()
                                                     : temp_index.get_chain(next_index).forward_loops.front();
                    if (loop_distance != std::numeric_limits<size_t>::max() &&
                        visited_nodes.count(make_pair(next_index, !next_rev)) == 0 &&
                        arriving_nid != temp_snarl_record.start_node_id &&
                        arriving_nid != temp_snarl_record.end_node_id
                        ) {
                        //If the next node can loop back on itself, then add the next node in the opposite direction
                        const auto& nchain = temp_index.get_chain(next_index);
                        nid_t boundary_id = next_rev ? nchain.end_node_id : nchain.start_node_id;
                        size_t boundary_len = temp_index.get_node(
                            {SnarlDistanceIndex::TEMP_NODE, boundary_id}).node_length;
                        size_t next_node_len = loop_distance + 2 * boundary_len;
                        queue.push(make_pair(SnarlDistanceIndex::sum(current_distance, next_node_len),
                                       make_pair(next_index, !next_rev)));
                    }
                }
#ifdef debug_distance_indexing
                cerr << "        reached child " << temp_index.structure_start_end_as_string(next_index) << " going "
                     << (next_rev ? "rev" : "fd") << " with distance " << current_distance << " for ranks " << start_rank << " " << next_rank << endl;
#endif
            });
        }
        if (is_internal_node && reachable_node_count != 1) {
            //If this is an internal node, then it must have only one edge for it to be a simple snarl
            temp_snarl_record.is_simple = false;
        }
    }

    /** Check the minimum length of the snarl passing through this node **/
    if (start_rank != 0 && start_rank != 1) {

        size_t child_max_length = start_index.first == SnarlDistanceIndex::TEMP_NODE 
            ? temp_index.get_node(start_index).node_length
            : temp_index.get_chain(start_index).max_length;
        //The distance through the whole snarl traversing this node forwards
        //(This might actually be traversing it backwards but it doesn't really matter)

        size_t dist_start_left = start_index.first == SnarlDistanceIndex::TEMP_NODE 
            ? temp_index.get_node(start_index).distance_left_start
            : temp_index.get_chain(start_index).distance_left_start;
        size_t dist_end_right = start_index.first == SnarlDistanceIndex::TEMP_NODE 
            ? temp_index.get_node(start_index).distance_right_end
            : temp_index.get_chain(start_index).distance_right_end;
        size_t dist_start_right =  start_index.first == SnarlDistanceIndex::TEMP_NODE 
            ? temp_index.get_node(start_index).distance_right_start
            : temp_index.get_chain(start_index).distance_right_start;
        size_t dist_end_left = start_index.first == SnarlDistanceIndex::TEMP_NODE 
            ? temp_index.get_node(start_index).distance_left_end
            : temp_index.get_chain(start_index).distance_left_end;

        size_t snarl_length_fd = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                dist_start_left, dist_end_right),child_max_length);
        //The same thing traversing this node backwards
        size_t snarl_length_rev = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                dist_start_right, dist_end_left), child_max_length);
        //The max that isn't infinite
        size_t max_length = 
            snarl_length_rev == std::numeric_limits<size_t>::max() 
            ? snarl_length_fd 
            : (snarl_length_fd == std::numeric_limits<size_t>::max() 
                    ? snarl_length_rev 
                    : std::max(snarl_length_rev, snarl_length_fd));
        if (max_length != std::numeric_limits<size_t>::max()) {
            temp_snarl_record.max_length = std::max(temp_snarl_record.max_length, max_length);
        }
        if ( temp_snarl_record.is_simple && 
            ! ((dist_start_left == 0 && dist_end_right == 0 && dist_end_left == std::numeric_limits<size_t>::max() && dist_start_right == std::numeric_limits<size_t>::max() ) || 
               (dist_start_left == std::numeric_limits<size_t>::max() && dist_end_right == std::numeric_limits<size_t>::max() && dist_end_left == 0 && dist_start_right == 0 ))){
            //If the snarl is simple, double check that this node is actually simple: that it can only be traversed going
            //across the nsarl
            temp_snarl_record.is_simple = false;
        }
    }
}

} // namespace vg
