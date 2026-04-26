//#define debug_subgraph

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {

void subgraph_in_distance_range(const SnarlDistanceIndex& distance_index, const Path& path, const HandleGraph* super_graph, size_t min_distance,
                                        size_t max_distance, std::unordered_set<nid_t>& subgraph, bool look_forward){

    //The position we're starting from - either the start or end of the path
    pos_t start_pos;
    size_t node_len;
    if (look_forward ){
        start_pos = initial_position(path);
        node_len = super_graph->get_length(super_graph->get_handle(get_id(start_pos)));
    } else {
        start_pos = final_position(path);
        node_len = super_graph->get_length(super_graph->get_handle(get_id(start_pos)));
        start_pos = reverse_base_pos(start_pos, node_len);
    }
    pair<nid_t, bool> traversal_start = std::make_pair(get_id(start_pos), get_is_rev(start_pos));

#ifdef debug_subgraph
cerr << endl << "Find subgraph in distance range " << min_distance << " to " << max_distance << endl;
cerr << "Start positon: "<< start_pos << endl;
#endif
    //The distance from the position to the ends of the current node(/snarl/chain)
    size_t current_distance_left = is_rev(start_pos) ? node_len - get_offset(start_pos) : std::numeric_limits<size_t>::max() ;
    size_t current_distance_right = is_rev(start_pos) ? std::numeric_limits<size_t>::max() : node_len - get_offset(start_pos) ;

    //Graph node of the start and end of the current node(/snarl/chain) pointing out
    net_handle_t current_net = distance_index.get_node_net_handle(get_id(start_pos));
    net_handle_t parent = distance_index.start_end_traversal_of(distance_index.get_parent(current_net));

    //The id and orientation of nodes that are too close and should be avoided
    hash_set<pair<id_t, bool>> seen_nodes;
    //Nodes that we want to start a search from - the distance is smaller or equal to than min_distance but
    //we can't walk out any further along the snarl tree without exceeding it
    //The distance is the distance from the start position to the beginning (or end if its backwards) of the node,
    //including the position
    vector<pair<handle_t, size_t>> search_start_nodes;

    if (((current_distance_left != std::numeric_limits<size_t>::max() && current_distance_left > min_distance) ||
           (current_distance_right != std::numeric_limits<size_t>::max() && current_distance_right > min_distance)) ||
         (distance_index.is_trivial_chain(parent) 
            && distance_index.distance_in_parent(distance_index.get_parent(parent), parent, distance_index.flip(parent)) == 0
            && node_len*2 > min_distance)) {
        //If the distance to either end of the node is within the range
        //Or of there is a loop on the node ( a duplication of just the node) and the node length would put one loop in the distance range

        //Add this node to the subgraph
        subgraph.emplace(get_id(start_pos));

        handle_t start = is_rev(start_pos) ? distance_index.get_handle(distance_index.flip(current_net), super_graph)
                                           : distance_index.get_handle(current_net, super_graph); 

        //Add any node one step out from this one to search_start_nodes
        super_graph->follow_edges(start, 
                false, [&](const handle_t& next_handle) {
            search_start_nodes.emplace_back(next_handle, is_rev(start_pos) ? current_distance_left : current_distance_right);
        });

        //Search for reachable nodes
        subgraph_in_distance_range_walk_graph(super_graph, min_distance, max_distance, subgraph, search_start_nodes, seen_nodes, traversal_start); 

        return;
    }


    
    while (!distance_index.is_root(parent)) {
#ifdef debug_subgraph
        cerr << "At child " << distance_index.net_handle_as_string(current_net) << " with distances " << current_distance_left << " " << current_distance_right << endl;
        cerr << "Parent is " << distance_index.net_handle_as_string(parent) << " at offset " << SnarlDistanceIndex::get_record_offset(parent) << endl;
#endif

        size_t max_parent_length = distance_index.maximum_length(parent);


        //Distances to get to the ends of the parent
        size_t distance_start_left = SnarlDistanceIndex::sum(current_distance_left,
                    distance_index.distance_to_parent_bound(parent, true, distance_index.flip(current_net)));
        size_t distance_start_right = SnarlDistanceIndex::sum(current_distance_right,
                     distance_index.distance_to_parent_bound(parent, true, current_net));
        size_t distance_end_left = SnarlDistanceIndex::sum(current_distance_left,
                    distance_index.distance_to_parent_bound(parent, false, distance_index.flip(current_net)));
        size_t distance_end_right = SnarlDistanceIndex::sum(current_distance_right,
                     distance_index.distance_to_parent_bound(parent, false, current_net));

        if ((current_distance_right != std::numeric_limits<size_t>::max() && current_distance_right >= min_distance)
            || (current_distance_left != std::numeric_limits<size_t>::max() && current_distance_left >= min_distance)
            || (distance_start_right != std::numeric_limits<size_t>::max() && distance_start_right>= min_distance)
            || (distance_end_right != std::numeric_limits<size_t>::max() && distance_end_right >= min_distance) 
            || (distance_start_left != std::numeric_limits<size_t>::max() && distance_start_left >= min_distance)
            || (distance_end_left != std::numeric_limits<size_t>::max() && distance_end_left >= min_distance)
            || (max_parent_length != std::numeric_limits<size_t>::max() && max_parent_length >= min_distance)) {
            //If the min distance will be exceeded within this parent, then start a search from the ends of this child

            if (distance_index.is_snarl(parent)) {
                //If this is the child of a snarl, then just traverse from the end of the node
#ifdef debug_subgraph
                cerr << "Start search in parent " << distance_index.net_handle_as_string(parent);
#endif
                if (current_distance_left != std::numeric_limits<size_t>::max() ){
                    //If we can go left
                    net_handle_t bound = distance_index.is_node(current_net) ? distance_index.flip(current_net) 
                                : distance_index.get_bound(current_net, false, false);
                    if (distance_index.is_sentinel(bound)) {
                        bound = distance_index.get_node_from_sentinel(bound);
                    }
                    handle_t current_node = distance_index.get_handle(bound, super_graph);
                    //Add everything immediately after the left bound of this node/chain
                    super_graph->follow_edges(distance_index.get_handle(bound, super_graph),
                            false, [&](const handle_t& next_handle) {
                        seen_nodes.erase(make_pair(super_graph->get_id(next_handle), super_graph->get_is_reverse(next_handle)));
                        search_start_nodes.emplace_back(next_handle,current_distance_left);

                    });

#ifdef debug_subgraph
                    cerr << " going left from " << super_graph->get_id(current_node) << (super_graph->get_is_reverse(current_node) ? "rev " : "fd ") ;
#endif
                } 
                if (current_distance_right != std::numeric_limits<size_t>::max()) {
                    //If we can go right
                    net_handle_t bound = distance_index.is_node(current_net) ? current_net 
                                : distance_index.get_bound(current_net, true, false);
                    if (distance_index.is_sentinel(bound)) {
                        bound = distance_index.get_node_from_sentinel(bound);
                    }
                    handle_t current_node = distance_index.get_handle(bound, super_graph);

                    //Add everything immediately after the right bound of this node/chain
                    super_graph->follow_edges(distance_index.get_handle(bound, super_graph),
                            false, [&](const handle_t& next_handle) {
                        seen_nodes.erase(make_pair(super_graph->get_id(next_handle),super_graph->get_is_reverse(next_handle)));
                        search_start_nodes.emplace_back(next_handle, current_distance_right);
                    });

#ifdef debug_subgraph
                    cerr << " going right from " << super_graph->get_id(current_node) << (super_graph->get_is_reverse(current_node) ? "rev " : "fd ");
#endif
                }
#ifdef debug_subgraph
                cerr << endl;
#endif
            } else {
#ifdef debug_subgraph
                cerr << "Start search along parent chain " << distance_index.net_handle_as_string(parent);
#endif
                //If this is the child of a chain, then traverse along the chain
                if (current_distance_left != std::numeric_limits<size_t>::max()) {
                    subgraph_in_distance_range_walk_across_chain (distance_index, super_graph,  subgraph, 
                        distance_index.flip(current_net), current_distance_left, search_start_nodes, seen_nodes, min_distance, max_distance, false);
                }
                if (current_distance_right != std::numeric_limits<size_t>::max()) {
                    subgraph_in_distance_range_walk_across_chain (distance_index, super_graph,  subgraph, 
                        current_net, current_distance_right, search_start_nodes, seen_nodes, min_distance, max_distance, false);
                }
            }
            subgraph_in_distance_range_walk_graph(super_graph, min_distance, max_distance, subgraph, search_start_nodes, seen_nodes, traversal_start); 
            return;
        } else if (distance_index.is_snarl(parent)){
#ifdef debug_subgraph
            cerr << "Parent is a snarl of handle type " << SnarlDistanceIndex::get_handle_type(parent) << " at offset " << SnarlDistanceIndex::get_record_offset(parent) << endl;
#endif
            //TODO: This might be overkill. It prevents us from adding nodes that shouldn't be in the subgraph, but might be too slow
            //If we don't check the other direction, go through the loop and add everything whose distance is lower than the minimum
            //to seen_nodes
            vector<pair<handle_t, size_t>> loop_handles_to_check;
            handle_t start_out = distance_index.get_handle(distance_index.get_bound(parent, false, false), super_graph);
            handle_t end_out = distance_index.get_handle(distance_index.get_bound(parent, true, false), super_graph);
            if (current_distance_left != std::numeric_limits<size_t>::max()) {
                loop_handles_to_check.emplace_back(distance_index.get_handle(distance_index.get_bound(current_net, false, false), super_graph), current_distance_left);
            }
            if (current_distance_right != std::numeric_limits<size_t>::max()) {
                loop_handles_to_check.emplace_back(distance_index.get_handle(distance_index.get_bound(current_net, true, false), super_graph), current_distance_right);
            }
            while (!loop_handles_to_check.empty()) {
                handle_t current_loop_handle = loop_handles_to_check.back().first;
                size_t current_loop_distance = loop_handles_to_check.back().second;
                loop_handles_to_check.pop_back();

                //Add to seen_nodes
                seen_nodes.emplace(super_graph->get_id(current_loop_handle), super_graph->get_is_reverse(current_loop_handle));

                //Walk one step out from this node
                super_graph->follow_edges(current_loop_handle, false, [&](const handle_t& next_handle) {
                    //If the next node is close enough and isn't exiting the snarl, then add it to stack
                    size_t new_distance = SnarlDistanceIndex::sum(current_loop_distance, super_graph->get_length(next_handle));
                    if (new_distance < min_distance && next_handle != start_out && next_handle != end_out) {
                        loop_handles_to_check.emplace_back(next_handle, new_distance);
                    }
                });
            }
        } else if (distance_index.is_chain(parent)) {
#ifdef debug_subgraph
            cerr << "Parent is a chain of handle type " << SnarlDistanceIndex::get_handle_type(parent) << " at offset " << SnarlDistanceIndex::get_record_offset(parent) << endl;
#endif
            //TODO: This is probably also overkill - walk a chain if there is a viable loop
            size_t distance_loop_right = distance_index.distance_in_parent(parent, current_net, current_net, super_graph, max_distance);
            size_t distance_loop_left =  distance_index.distance_in_parent(parent, distance_index.flip(current_net), distance_index.flip(current_net), super_graph, max_distance);
            if ((current_distance_left != std::numeric_limits<size_t>::max() && distance_loop_left != std::numeric_limits<size_t>::max()) ||
                (current_distance_right != std::numeric_limits<size_t>::max() && distance_loop_right != std::numeric_limits<size_t>::max())) {
                //If there is a loop that we can take, then take it
                if (current_distance_left != std::numeric_limits<size_t>::max()) {
                    subgraph_in_distance_range_walk_across_chain (distance_index, super_graph,  subgraph, 
                        distance_index.flip(current_net), current_distance_left, search_start_nodes, seen_nodes, min_distance, max_distance, false);
                }
                if (current_distance_right != std::numeric_limits<size_t>::max()) {
                    subgraph_in_distance_range_walk_across_chain (distance_index, super_graph,  subgraph, 
                        current_net, current_distance_right, search_start_nodes, seen_nodes, min_distance, max_distance, false);
                }
                subgraph_in_distance_range_walk_graph(super_graph, min_distance, max_distance, subgraph, search_start_nodes, seen_nodes, traversal_start); 
                return;
            }
        }

        //Remember the bounds of this child so we don't return to it
        if (current_distance_left != std::numeric_limits<size_t>::max() ){
            //If we can go left
            net_handle_t bound = distance_index.is_node(current_net) ? distance_index.flip(current_net) 
                        : distance_index.get_bound(current_net, false, false);
            if (distance_index.is_sentinel(bound)) {
                bound = distance_index.get_node_from_sentinel(bound);
            }
            handle_t current_node = distance_index.get_handle(bound, super_graph);
            seen_nodes.emplace(super_graph->get_id(current_node), super_graph->get_is_reverse(current_node));
        }
        if (current_distance_right != std::numeric_limits<size_t>::max()) {
            //If we can go right
            net_handle_t bound = distance_index.is_node(current_net) ? current_net 
                        : distance_index.get_bound(current_net, true, false);
            if (distance_index.is_sentinel(bound)) {
                bound = distance_index.get_node_from_sentinel(bound);
            }
            handle_t current_node = distance_index.get_handle(bound, super_graph);
            seen_nodes.emplace(super_graph->get_id(current_node), super_graph->get_is_reverse(current_node));
        }

        current_distance_left = std::min(distance_start_left, distance_start_right);
        current_distance_right = std::min(distance_end_left, distance_end_right);

        current_net = std::move(parent);
        parent = distance_index.canonical(distance_index.get_parent(current_net));
    }
    if (current_distance_left <= min_distance) {
#ifdef debug_subgraph
        cerr << "Adding the end of a child of the root " << distance_index.net_handle_as_string(distance_index.get_bound(current_net, false, false)) << " with distance " << current_distance_left << endl;
#endif

        handle_t bound = distance_index.get_handle(distance_index.get_bound(current_net, false, false), super_graph);
        search_start_nodes.emplace_back(bound, current_distance_left);
    }
    if (current_distance_right <= min_distance) {
#ifdef debug_subgraph
        cerr << "Adding the end of a child of the root " << distance_index.net_handle_as_string(distance_index.get_bound(current_net, false, false)) << " with distance " << current_distance_right << endl;
#endif
        handle_t bound = distance_index.get_handle(distance_index.get_bound(current_net, true, false), super_graph);
        search_start_nodes.emplace_back(bound,current_distance_right);
    }
    subgraph_in_distance_range_walk_graph(super_graph, min_distance, max_distance, subgraph, search_start_nodes, seen_nodes, traversal_start); 

    return;
}


///Helper for subgraph_in_distance_range
///Given starting handles in the super graph and the distances to each handle (including the start position and
//the first position in the handle), add all nodes within the distance range, excluding nodes in seen_nodes
void subgraph_in_distance_range_walk_graph(const HandleGraph* super_graph, size_t min_distance, size_t max_distance,
                        std::unordered_set<nid_t>& subgraph, vector<pair<handle_t, size_t>>& start_nodes,
                        hash_set<pair<nid_t, bool>>& seen_nodes, const pair<nid_t, bool>& traversal_start) {
#ifdef debug_subgraph
    cerr << "Starting search from nodes " << endl;
    for (auto& start_handle : start_nodes) {
        cerr << "\t" << super_graph->get_id(start_handle.first) << " " << super_graph->get_is_reverse(start_handle.first)
             << " with distance " << start_handle.second << endl;
    }
#endif

    //Order based on the distance to the position (handle)
    auto cmp =  [] (const pair<handle_t, size_t> a, const pair<handle_t, size_t> b ) {
            return a.second > b.second;
        };
    priority_queue< pair<handle_t, size_t>, vector<pair<handle_t, size_t>>, decltype(cmp)> next_handles (cmp);
    for (auto& start_handle : start_nodes) {
        next_handles.emplace(start_handle);
    }
    bool first_node = true;

    while (next_handles.size() > 0) {
        //Traverse the graph, adding nodes if they are within the range
        handle_t curr_handle=next_handles.top().first;
        size_t curr_distance=next_handles.top().second;
        next_handles.pop();
#ifdef debug_subgraph
        cerr << "At node " << super_graph->get_id(curr_handle) << " " << super_graph->get_is_reverse(curr_handle) << " with distance " << curr_distance << endl;
#endif
        if (seen_nodes.count(make_pair(super_graph->get_id(curr_handle), super_graph->get_is_reverse(curr_handle))) == 0) {
            seen_nodes.emplace(super_graph->get_id(curr_handle), super_graph->get_is_reverse(curr_handle));

            size_t node_len = super_graph->get_length(curr_handle);
            size_t curr_distance_end = SnarlDistanceIndex::sum(curr_distance, node_len)-1;
            if ((curr_distance >= min_distance && curr_distance <= max_distance) ||
                 (curr_distance_end >= min_distance && curr_distance_end <= max_distance) ||
                 (curr_distance <= min_distance && curr_distance_end >= max_distance)) {
#ifdef debug_subgraph
                cerr << "\tadding node " << super_graph->get_id(curr_handle) << " " << super_graph->get_is_reverse(curr_handle) << " with distance "
                     << curr_distance << " and node length " << node_len << endl;
#endif
                subgraph.insert(super_graph->get_id(curr_handle));

            }
#ifdef debug_subgraph
            else {
                cerr << "\tdisregarding node " << super_graph->get_id(curr_handle) << " " << super_graph->get_is_reverse(curr_handle)
                     << " with distance " << curr_distance << " and node length " << node_len << endl;
            }
#endif
            curr_distance = SnarlDistanceIndex::sum(node_len, curr_distance);

            //If the end of this node is still within the range, add the next nodes that are within
            //Also check that the node we're currently at isn't the start node
            if (SnarlDistanceIndex::minus(curr_distance,1) <= max_distance) {
                super_graph->follow_edges(curr_handle, false, [&](const handle_t& next) {
                    nid_t next_id = super_graph->get_id(next);
                    if (seen_nodes.count(make_pair(next_id, super_graph->get_is_reverse(next))) == 0) {
                        next_handles.emplace(next, curr_distance);
                    } 
                    return true;
                });
            }
            first_node = false;
        } 
#ifdef debug_subgraph 
        else {
            cerr << "\tthe node was already seen" << endl;
        }
#endif

    }

#ifdef debug_subgraph
    cerr << "Subgraph has nodes: ";
    for (const nid_t& node : subgraph) {
        cerr << node << ", ";
    }
    cerr << endl;
#endif
    return;
}
//helper function to walk along a chain from the current node until the distance traversed
//exceeds the minimum limit. Add the node just before this happens to search_start_nodes
void subgraph_in_distance_range_walk_across_chain (const SnarlDistanceIndex& distance_index, const HandleGraph* super_graph,
        std::unordered_set<nid_t>& subgraph, net_handle_t current_node, 
        size_t current_distance, vector<pair<handle_t, size_t>>& search_start_nodes, hash_set<pair<nid_t, bool>>& seen_nodes, 
        const size_t& min_distance, const size_t& max_distance, bool checked_loop){
#ifdef debug_subgraph
    cerr << "Walk along parent chain " << distance_index.net_handle_as_string(distance_index.get_parent(current_node)) << " from " << distance_index.net_handle_as_string(current_node) << " with " << current_distance << endl;
#endif
    if (distance_index.is_trivial_chain(distance_index.get_parent(current_node))){
        return;
    }
    bool finished_chain = false;
    bool added_nodes = false; //Did we start a search? if not, add the last node in the chain
    while (current_distance <= min_distance && !finished_chain) {
        finished_chain = distance_index.follow_net_edges(current_node, super_graph, false, 
            [&](const net_handle_t& next) {
                size_t next_length = distance_index.minimum_length(next);
                //If the next child is a snarl, then the distance to loop in the snarl
                if (distance_index.is_snarl(next)) {
                    net_handle_t bound_fd = distance_index.get_bound(next, distance_index.ends_at(next) == SnarlDistanceIndex::START, true);
                    size_t next_loop = distance_index.distance_in_parent(next, bound_fd, bound_fd, super_graph, max_distance);
                    if (!checked_loop && next_loop != std::numeric_limits<size_t>::max()) {
#ifdef debug_subgraph
                        cerr << "\tsnarl loops so also check the other direction" << endl;
#endif
                        //If we haven't yet checked the chain in the other direction and this snarl allows us to loop
                        if ( SnarlDistanceIndex::sum(next_loop, current_distance) != std::numeric_limits<size_t>::max()  &&
                             SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(next_loop, 
                                                                             current_distance), 
                                                                             distance_index.node_length(current_node)) >= min_distance) {
#ifdef debug_subgraph
                            cerr << "\t\t add the current node" << endl;
#endif
                            //If the loop will put us over the edge, then start from the current node
                            super_graph->follow_edges(distance_index.get_handle(current_node, super_graph), false, [&](const handle_t& next_handle) {
                                search_start_nodes.emplace_back(next_handle,current_distance);
                            });
                            return true;
                        } else {
                            //Otherwise, switch direction in the chain and walk along it again
                            subgraph_in_distance_range_walk_across_chain(distance_index, super_graph, subgraph, distance_index.flip(current_node),
                                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(current_distance, 
                                                                                    next_loop), 
                                                                                    distance_index.node_length(current_node)), 
                                    search_start_nodes, seen_nodes, min_distance, max_distance, true);
                            checked_loop = true;
                        }
                    }
                    if (next_loop != std::numeric_limits<size_t>::max()){
                        //TODO: This might be overkill. It prevents us from adding nodes that shouldn't be in the subgraph, but might be too slow
                        //If we don't check the other direction, go through the loop and add everything whose distance is lower than the minimum
                        //to seen_nodes
                        vector<pair<handle_t, size_t>> loop_handles_to_check;
                        handle_t start_out = distance_index.get_handle(distance_index.get_bound(next, false, false), super_graph);
                        handle_t end_out = distance_index.get_handle(distance_index.get_bound(next, true, false), super_graph);
                        loop_handles_to_check.emplace_back(distance_index.get_handle(bound_fd, super_graph), current_distance);
                        while (!loop_handles_to_check.empty()) {
                            handle_t current_loop_handle = loop_handles_to_check.back().first;
                            size_t current_loop_distance = loop_handles_to_check.back().second;
                            loop_handles_to_check.pop_back();

                            //Add to seen_nodes
                            seen_nodes.emplace(super_graph->get_id(current_loop_handle), super_graph->get_is_reverse(current_loop_handle));

                            //Walk one step out from this node
                            super_graph->follow_edges(current_loop_handle, false, [&](const handle_t& next_handle) {
                                //If the next node is close enough and isn't exiting the snarl, then add it to stack
                                size_t new_distance = SnarlDistanceIndex::sum(current_loop_distance, super_graph->get_length(next_handle));
                                if (new_distance < min_distance && next_handle != start_out && next_handle != end_out) {
                                    loop_handles_to_check.emplace_back(next_handle, new_distance);
                                }
                            });
                        }

                    }
                }
                size_t next_max_length = distance_index.maximum_length(next);
#ifdef debug_subgraph
                cerr << "\tnext node: " << distance_index.net_handle_as_string(next) << " with distance " << current_distance << " and min and max lengths " << next_length << " " << next_max_length << endl;
#endif
                if (( SnarlDistanceIndex::sum(next_max_length, current_distance) != std::numeric_limits<size_t>::max()  &&
                     SnarlDistanceIndex::sum(next_max_length, current_distance) >= min_distance)){
                    if (distance_index.is_node(next)) {
                        size_t curr_distance_end = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(next_max_length, current_distance),1);
                        //If its a node that puts us over, add the node to the subgraph, then start the search from that node
#ifdef debug_subgraph
                        cerr << "\t\tAdding node from a chain " << distance_index.net_handle_as_string(next) << " with distance " << current_distance << endl;
#endif
                        if ((current_distance >= min_distance && current_distance <= max_distance) ||
                             (curr_distance_end >= min_distance && curr_distance_end <= max_distance) ||
                             (current_distance <= min_distance && curr_distance_end >= max_distance)) {
                            subgraph.emplace(distance_index.node_id(next));
                        }
                        super_graph->follow_edges(distance_index.get_handle(next, super_graph), false, [&](const handle_t& next_handle) {
                            search_start_nodes.emplace_back(next_handle, SnarlDistanceIndex::sum(current_distance, next_length));
                            seen_nodes.erase(make_pair(super_graph->get_id(next_handle), super_graph->get_is_reverse(next_handle)));
                        });
                    } else {
                        //If it's a snarl, then we'll start from the last node
#ifdef debug_subgraph
                        cerr << "\t\tAdding node from a chain " << distance_index.net_handle_as_string(next) << " with distance " << current_distance << endl;
#endif
                        super_graph->follow_edges(distance_index.get_handle(current_node, super_graph), false, [&](const handle_t& next_handle) {
                            search_start_nodes.emplace_back(next_handle,current_distance);
                            seen_nodes.erase(make_pair(super_graph->get_id(next_handle), super_graph->get_is_reverse(next_handle)));
                        });
                    }
                    //If we added something, stop traversing the chain
                    added_nodes = true;
                    return true;
                } else if (distance_index.is_node(next)) {
                    seen_nodes.emplace(distance_index.node_id(next), distance_index.ends_at(next) == SnarlDistanceIndex::START);
                }
                current_node = next;
                current_distance = SnarlDistanceIndex::sum(next_length, current_distance);
                if (current_distance > max_distance) {
                    added_nodes = true;
                    return true;
                } else {
                    return false;
                }
        }); 
    }
    if (!added_nodes && current_distance <= max_distance) {
        //If we haven't added anything and haven't exceeded the distance limit, then start from the end of the chain
        handle_t bound = distance_index.get_handle(current_node, super_graph);

        super_graph->follow_edges(bound, false, [&](const handle_t& next_handle) {
            search_start_nodes.emplace_back(next_handle,current_distance);
            seen_nodes.erase(make_pair(super_graph->get_id(next_handle), super_graph->get_is_reverse(next_handle)));
        });
        //seen_nodes.erase(make_pair(super_graph->get_id(bound), super_graph->get_is_reverse(bound)));
        //search_start_nodes.emplace_back( bound, current_distance);
    }
};

} // namespace vg
