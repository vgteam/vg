#include <unordered_set>
#include "cactus.hpp"
#include "vg.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/weakly_connected_components.hpp"
#include "algorithms/find_shortest_paths.hpp"

extern "C" {
#include "sonLib.h"
#include "stCactusGraphs.h"
}

namespace vg {

using namespace std;

// Cactus


// Step 1) Get undirected adjacency connected components of VG *sides*
// (note we can't put NodeSides in header, so leave this function static)
// Here and elsewhere we adopt the convention that a forward handle indicates
// the end of a node
typedef unordered_set<handle_t> HandleSet;
typedef unordered_map<handle_t, int> Handle2Component;
static void compute_side_components(const HandleGraph& graph,
                                    vector<HandleSet>& components,
                                    Handle2Component& handle_to_component) {
    
    // if a handle's side hasn't been added to a component yet, traverses it's component
    // and adds it to the return structures
    function<void(const handle_t& handle)> add_handle = [&](const handle_t& handle) {
        // have we already visited this side?
        if (!handle_to_component.count(handle)) {
            // make a new component
            int component_id = components.size();
            components.emplace_back();
            auto& component = components.back();
            
            // mark the current node as being on this component
            handle_to_component[handle] = component_id;
            component.insert(handle);
            
            // set up structures to do a BFS traversal
            queue<handle_t> q;
            q.push(handle);
            
            while (!q.empty()) {
                handle_t next_handle = q.front();
                q.pop();
                
                graph.follow_edges(next_handle, false, [&](const handle_t& traversed_handle) {
                    // to get the handle corresponding to the side we're touching, we have to
                    // flip the traversal around
                    handle_t adjacent_side = graph.flip(traversed_handle);
                    
                    // if necessary, enqueue and assign the side to a component
                    if (!handle_to_component.count(adjacent_side)) {
                        handle_to_component[adjacent_side] = component_id;
                        component.insert(adjacent_side);
                        q.push(adjacent_side);
                    }
                    return true;
                });
            }
        }
    };

    graph.for_each_handle([&](const handle_t& handle) {
        // component for end of node
        add_handle(handle);
        // component for start of node
        add_handle(graph.flip(handle));
    });
    
}

void* mergeNodeObjects(void* a, void* b) {
    // One of the objects is going to get returned, and the other one is going
    // to get freed (since the graph is supposed to own them all).
    id_t* to_return;
    id_t* to_free;
    
    if(*(id_t*)a < *(id_t*)b) {
        to_return = (id_t*)a;
        to_free = (id_t*)b;
    } else {
        to_free = (id_t*)a;
        to_return = (id_t*)b;
    }
    
#ifdef debug
    cerr << "Free object " << to_free << " = " << *to_free << " and keep object " << to_return << " = " << *to_return << endl;
#endif
    
    // Free the thing we aren't keeping
    free(to_free);
    
    // Return the one we are
    return (void*)to_return;
}

/**
 * Get the bridge ends that form boundary pairs with edgeEnd1, using the given
 * getBridgeEdgeEndsToBridgeNodes hash map. Duplicated from the pinchesAndCacti
 * tests.
 */
void getReachableBridges2(stCactusEdgeEnd *edgeEnd1,
        stHash *bridgeEndsToBridgeNodes, stList *bridgeEnds) {
    assert(edgeEnd1->link == NULL); // is a bridge

    // The node in the bridge graph incident with edgeEnd1
    stBridgeNode *bNode = (stBridgeNode*)stHash_search(bridgeEndsToBridgeNodes, edgeEnd1);

    // Walk from bNode to all the reachable nodes
    stSetIterator *endIt = stSet_getIterator(bNode->bridgeEnds);
    stCactusEdgeEnd *edgeEnd2;
    while((edgeEnd2 = (stCactusEdgeEnd*)stSet_getNext(endIt)) != NULL) {
        if(edgeEnd2 != edgeEnd1) {
            stList_append(bridgeEnds, edgeEnd2);
            getReachableBridges2(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2), bridgeEndsToBridgeNodes, bridgeEnds);
        }
    }
    stSet_destructIterator(endIt);
}

/**
 * Get the bridge ends that form boundary pairs with edgeEnd1.
 * Duplicated from the pinchesAndCacti tests.
 */
void getReachableBridges(stCactusEdgeEnd *edgeEnd1, stList *bridgeEnds) {

    // Get bridge graph and map of bridge ends to bridge nodes in the bridge graph
    stBridgeGraph *bGraph = stBridgeGraph_getBridgeGraph(edgeEnd1->node);
    stHash *bridgeEndsToBridgeNodes = stBridgeGraph_getBridgeEdgeEndsToBridgeNodesHash(bGraph);

    // Do DFS to get reachable bridges
    getReachableBridges2(edgeEnd1, bridgeEndsToBridgeNodes, bridgeEnds);

    // Cleanup
    stHash_destruct(bridgeEndsToBridgeNodes);
    stBridgeGraph_destruct(bGraph);
}

/**
 * Finds an arbitrary pair of telomeres in a Cactus graph, which are are either
 * a pair of bridge edge ends or a pair of chain edge ends, oriented such that
 * they form a pair of boundaries.
 *
 * Mostly copied from the pinchesAndCacti unit tests.
 */
void addArbitraryTelomerePair(vector<stCactusEdgeEnd*> ends, stList *telomeres) {

    // If empty graph, print warning and exit
    if(ends.empty()) {
        throw runtime_error("Empty graph, no telomeres to select");
    }

    // Pick an arbitrary edge end
    stCactusEdgeEnd *edgeEnd1 = ends.front();

    // Now find a compatible end
    stCactusEdgeEnd *edgeEnd2;

    // If a chain end
    if(edgeEnd1->link != NULL) {
        // Get another elligible edge end in the chain that forms a pair of boundaries
        edgeEnd2 = edgeEnd1->link;
        
        // TODO: we could also follow edgeEnd2->otherEdgeEnd->link and get another
        // valid option, until we start to circle around and repeat.
    } else {
        // Else, is a bridge end.
        // Get the other bridges in the subtree.
        stList *bridgeEnds = stList_construct();
        getReachableBridges(edgeEnd1, bridgeEnds);
        stList_append(bridgeEnds, edgeEnd1); // Case it is a unary top-level snarl

        // Get an arbitrary list member
        edgeEnd2 = (stCactusEdgeEnd*)stList_get(bridgeEnds, 0);

        // Cleanup
        stList_destruct(bridgeEnds);
    }

    // Add to telomeres
    stList_append(telomeres, edgeEnd1);
    stList_append(telomeres, edgeEnd2);
}

// Step 2) Make a Cactus Graph. Returns the graph and a list of paired
// cactusEdgeEnd telomeres, one after the other. Both members of the return
// value must be destroyed.
pair<stCactusGraph*, stList*> handle_graph_to_cactus(PathHandleGraph& graph, const unordered_set<string>& hint_paths) {

    // in a cactus graph, every node is an adjacency component.
    // every edge is a *vg* node connecting the component

    // start by identifying adjacency components
    vector<HandleSet> components;
    Handle2Component handle_to_component;
    compute_side_components(graph, components, handle_to_component);

    // map cactus nodes back to components
    vector<stCactusNode*> cactus_nodes(components.size());
    
    // create cactus graph
    stCactusGraph* cactus_graph = stCactusGraph_construct2(free, free);
    
    // copy each component into cactus node
    for (int i = 0; i < components.size(); ++i) {
        id_t* cactus_node_id = (id_t*)malloc(sizeof(id_t));
        *cactus_node_id = i;
        cactus_nodes[i] = stCactusNode_construct(cactus_graph, cactus_node_id);
        
#ifdef debug
        cerr << "Created cactus node " << cactus_nodes[i] << " with object " << cactus_node_id
            << " for component " << i << " with size " << components[i].size() << endl;
#endif
    }

    // Make edge for each vg node connecting two adjacency components. We also
    // keep track of the main cactusEdgeEnd we get for each node in its local
    // forward orientation.
    unordered_map<id_t, stCactusEdgeEnd*> edge_ends;
    for (int i = 0; i < components.size(); ++i) {
        // For each adjacency component
        for (auto side : components[i]) {
            // For every side in it
            
            // Work out the other side of that node
            handle_t other_side = graph.flip(side);
            // And what component it is in
            int j = handle_to_component[other_side];
            
            id_t node_id = graph.get_id(side);
            id_t other_node_id = graph.get_id(other_side);
            
            // by our convention, the forward strand traversal is the end
            bool is_end = !graph.get_is_reverse(side);
            bool other_is_end = !graph.get_is_reverse(other_side);
            
            if (!edge_ends.count(node_id)) {
                // If we haven't made the Cactus edge for this graph node yet
                
                // afraid to try to get C++ NodeSide class into C, so we copy to
                // equivalent struct
                CactusSide* cac_side1 = (CactusSide*)malloc(sizeof(CactusSide));
                CactusSide* cac_side2 = (CactusSide*)malloc(sizeof(CactusSide));
                cac_side1->node = node_id;
                cac_side1->is_end = is_end;
                cac_side2->node = other_node_id;
                cac_side2->is_end = other_is_end;
#ifdef debug
                cerr << "Creating cactus edge for sides " << pb2json(to_visit(side)) << " -- " << pb2json(to_visit(other_side)) << ": "
                << i << " -> " << j << endl;
#endif
                
                // We get the cactusEdgeEnd corresponding to the side stored in side.
                // This may be either the left (if that NodeSide is a start), or the right (if that NodeSide is an end).
                stCactusEdgeEnd* cactus_edge = stCactusEdgeEnd_construct(cactus_graph,
                                                                         cactus_nodes[i],
                                                                         cactus_nodes[j],
                                                                         cac_side1,
                                                                         cac_side2);
                // Save the cactusEdgeEnd for the left side of the node.
                edge_ends[node_id] = is_end ? stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge) : cactus_edge;
            }
        }
    }

    // collapse 3-edge connected components to make actual cactus graph. 
    stCactusGraph_collapseToCactus(
        cactus_graph, mergeNodeObjects, cactus_nodes[0]);
        
    // Define a list of telomeres
    stList *telomeres = stList_construct();
    
    // Now we decide on telomere pairs.
    // We need one for each weakly connected component in the graph, so first we break into connected components.
    vector<unordered_set<id_t>> weak_components = algorithms::weakly_connected_components(&graph);
       
    // We also want a map so we can efficiently find which component a node lives in.
    unordered_map<id_t, size_t> node_to_component;
    for (size_t i = 0; i < weak_components.size(); i++) {
        if (weak_components[i].size() == 1) {
            // If we feed this through to Cactus it will crash.
            throw runtime_error("Cactus does not currently support finding snarls in a single-node connected component");
        }
    
        for (auto& id : weak_components[i]) {
            node_to_component[id] = i;
        }
    }
       
    // Then we find the heads and tails
    auto all_heads = algorithms::head_nodes(&graph);
    auto all_tails = algorithms::tail_nodes(&graph);
    
    // Alot them to components. We store tips in an inward-facing direction
    vector<unordered_set<handle_t>> component_tips(weak_components.size());
    for (auto& head : all_heads) {
        component_tips[node_to_component[graph.get_id(head)]].insert(head);
#ifdef debug
        cerr << "Found head " << graph.get_id(head) << " in component " << node_to_component[graph.get_id(head)] << endl;
#endif
    }
    for (auto& tail : all_tails) {
        component_tips[node_to_component[graph.get_id(tail)]].insert(graph.flip(tail));
#ifdef debug
        cerr << "Found tail " << graph.get_id(tail) << " in component " << node_to_component[graph.get_id(tail)] << endl;
#endif
    }
    
    // Assign path names to components
    vector<vector<string>> component_paths(weak_components.size());
    // Also get the path length.
    unordered_map<string, size_t> path_length;
    
    
    
    graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        
        if (graph.get_occurrence_count(path_handle) == 0) {
            // Not a real useful path, so skip it. Some alt paths used for
            // haplotype generation are empty.
            return;
        }
        
        string name = graph.get_path_name(path_handle);
        
        occurrence_handle_t occurrence_handle = graph.get_first_occurrence(occurrence_handle);
        
        auto component = node_to_component[graph.get_node_id(graph.get_occurrence(occurrence_handle))];
        component_paths[component].push_back(name);
        
        
#ifdef debug
        cerr << "Path " << name << " belongs to component " << component << endl;
#endif
        
        auto process_occurrence = [&](const occurrence_handle_t& occurrence_handle) {
            handle_t handle = graph.get_occurrence(occurrence_handle);
            path_length[name] += graph.get_length(handle);
            
            if (node_to_component[graph.get_node_id(handle)] != component) {
                // If we use a path like this to pick telomeres we will segfault Cactus.
                throw runtime_error("Path " + name + " spans multiple connected components!");
            }
        }
        
        while (graph.has_next_occurrence(occurrence_handle)) {
            process_occurrence(occurrence_handle);
            occurrence_handle = graph.get_next_occurrence(occurrence_handle);
        }
        process_occurrence(occurrence_handle);
        
#ifdef debug
        cerr << "\tPath " << name << " has length " << path_length[name] << endl;
#endif
    });
    
    // We'll also need the strongly connected components, in case the graph is cyclic.
    // This holds all the strongly connected components that live in each weakly connected component.
    vector<vector<set<id_t>>> component_strong_components(weak_components.size());
    for (auto& strong_component : graph.strongly_connected_components()) {
        // For each strongly connected component
        assert(!strong_component.empty());
        // Assign it to the weak comnponent that some node in it belongs to
        component_strong_components[node_to_component[*strong_component.begin()]].push_back(strong_component);
    }
    
    
    // OK, now we need to fill in the telomeres list with two telomeres per
    // component.
    
    // This function adds as telomeres a pair of inward-facing handles.
    auto add_telomeres = [&](const handle_t& left, const handle_t& right) {
#ifdef debug
        cerr << "Selected " << graph.get_id(left) << " " << graph.get_is_reverse(left) << " and "
            << graph.get_id(right) << " " << graph.get_is_reverse(right) << " as tips" << endl;
#endif
        
        stCactusEdgeEnd* end1 = edge_ends[graph.get_id(left)];
        // We need to add the interior side of the node, and our handle is reading inwards.
        // If we're reverse, add the node's local left. Otherwise, add the node's local right.
        stList_append(telomeres, graph.get_is_reverse(left) ? end1 : stCactusEdgeEnd_getOtherEdgeEnd(end1));
        stCactusEdgeEnd* end2 = edge_ends[graph.get_id(right)];
        stList_append(telomeres, graph.get_is_reverse(right) ? end2 : stCactusEdgeEnd_getOtherEdgeEnd(end2));
    };
    
    for (size_t i = 0; i < weak_components.size(); i++) {
        // For each weakly connected component
        
        // First priority is longest path in the hint set.
        // Next priority is two tips at the ends of the longest path that has
        // two tips at its ends.
        {
        
            string longest_path;
            // This is going to hold inward-facing handles to the tips we find.
            pair<handle_t, handle_t> longest_path_tips;
        
            for (bool require_hint : {true, false}) {
                // First try for a path that's in the hint set, then try for a path that might not be.
        
                size_t longest_path_length = 0;
#ifdef debug
                cerr << "Consider " << component_paths[i].size() << " paths for component " << i << endl;
#endif
                for (auto& path_name : component_paths[i]) {
                    // Look at each path
                    
                    if (require_hint && !hint_paths.count(path_name)) {
                        // Skip this one because it's not hinted
                        continue;
                    }
                    
                    auto& path_mappings = graph.paths.get_path(path_name);
                    
#ifdef debug
                    cerr << "\tPath " << path_name << " has " << path_mappings.size() << " mappings" << endl;
#endif
                    
                    // See if I can get two tips on its ends.
                    // Get the inward-facing start and end handles.
                    handle_t path_start = graph.get_handle(path_mappings.front().node_id(),
                        path_mappings.front().is_reverse());
                    handle_t path_end = graph.get_handle(path_mappings.back().node_id(),
                        !path_mappings.back().is_reverse());
                    
                    if (component_tips[i].count(path_start) && component_tips[i].count(path_end)) {
                        // This path ends in two tips so we can consider it
                        
                        if (path_length[path_name] > longest_path_length) {
                            // This is our new longest path between tips.
                            longest_path = path_name;
                            longest_path_tips = make_pair(path_start, path_end);
                            longest_path_length = path_length[path_name];
#ifdef debug
                            cerr << "\t\tNew longest path!" << endl;
#endif
                        } else {
#ifdef debug
                            cerr << "\t\tPath length of " << path_length[path_name] << " not longer than path "
                                << longest_path << " with " << longest_path_length << endl;
#endif
                        } 
                        
                    } else {
#ifdef debug
                        cerr << "\t\tPath " << path_name << " does not start and end with tips" << endl;
#endif
                    }
                }
                
                if (!longest_path.empty()) {
                    // We found something. Don't try again; we might be doing the hint pass and we don't want to clobber it.
                    break;
                }
                
            }
            
            if (!longest_path.empty()) {
#ifdef debug
                cerr << "Longest tip path is " << longest_path << endl;
#endif
                // We found something!
                add_telomeres(longest_path_tips.first, longest_path_tips.second);
                // Work on the next component.
                continue;
            } else {
#ifdef debug
                cerr << "No longest named path found between tips" << endl;
#endif
            }
            
        }
        
        // Otherwise, compute tip reachability and distance by Dijkstra's
        // Algorithm and pick the pair of reachable tips with the longest shortest
        // paths
        
        {
            // Track the distances between pairs of tips (or numeric_limits<size_t>::max() for unreachable).
            // Pairs are stored lowest-node-and-orientation first. Really we want a hashable unordered_pair...
            unordered_map<pair<handle_t, handle_t>, size_t> tip_distances;
            
            // Define a function to look up in the memo, or do Dijkstra if necessary
            // Takes a pair that must be sorted already.
            auto get_or_compute_distance = [&](pair<handle_t, handle_t>& key) {
                if (!tip_distances.count(key)) {
                    // If we don't know the distance, do a search out from one handle
                    
#ifdef debug
                    cerr << "Do Dijkstra traversal out from "
                        << graph.get_id(key.first) << " " << graph.get_is_reverse(key.first) << endl;
#endif
                    
                    unordered_map<handle_t, size_t> distances = algorithms::find_shortest_paths(&graph, key.first);
                    
                    for (auto& other_tip : component_tips[i]) {
                        // And save the distances for everything reachable or unreachable.
                        // This minimizes the number of searches we need to do.
                        
                        pair<handle_t, handle_t> other_key = make_pair(key.first, other_tip);
                        
                        if (graph.get_id(other_key.second) < graph.get_id(other_key.first) ||
                            (graph.get_id(other_key.second) == graph.get_id(other_key.first) &&
                            graph.get_is_reverse(other_key.second) < graph.get_is_reverse(other_key.first))) {
                            
                            // We need to put these tips the other way around
                            std::swap(other_key.first, other_key.second);
                        }
                                
                        if (distances.count(graph.flip(other_tip))) {
                            // Can get from tip 1 to this tip.
                            // We need to flip the tip because Dijkstra will be reading out of the graph and into the tip.
                            tip_distances[other_key] = distances[graph.flip(other_tip)];
                        } else {
                            // This tip is completely unreachable from tip 1
                            tip_distances[other_key] = numeric_limits<size_t>::max();
                        }
                    }
                }
                
                // Now either we can reach this tip or we can't. Either way it's in the memo.
                return tip_distances.at(key);
            }; 
            
            // Track the best pair of handles
            pair<handle_t, handle_t> furthest;
            // And their distance
            size_t furthest_distance = numeric_limits<size_t>::max();

            for (auto& tip1 : component_tips[i]) {
                for (auto& tip2 : component_tips[i]) {
                    // For each pair of tips
                    if (tip1 == tip2) {
                        // Skip unless they're distinct
                        continue;
                    }
           
                    // Make a pair to be the key
                    pair<handle_t, handle_t> key = make_pair(tip1, tip2);
                    if (graph.get_id(key.second) < graph.get_id(key.first) ||
                        (graph.get_id(key.second) == graph.get_id(key.first) &&
                        graph.get_is_reverse(key.second) < graph.get_is_reverse(key.first))) {
                        
                        // We need to put these tips the other way around
                        std::swap(key.first, key.second);
                    }
           
                    // Get or calculate their distance
                    auto tip_distance = get_or_compute_distance(key);
                    
                    if (tip_distance != numeric_limits<size_t>::max() &&
                        (tip_distance > furthest_distance || furthest_distance == numeric_limits<size_t>::max())) {
                        // If it's a new longer distance (but not infinite), between distinct tips, keep it
                        furthest_distance = tip_distance;
                        furthest = key;
                    }
                    
                }
            }
            
            if (furthest_distance == numeric_limits<size_t>::max()) {
                // We couldn't find anything with two different tips
#ifdef debug
                cerr << "No Dijkstra path found between distinct tips. Consider hairpins of one tip." << endl;
#endif

                // So now we will check for the same tip twice.
                for (auto& tip : component_tips[i]) {
                    // Make a key for two of each tip.
                    // No need to flip the symmetrical pair.
                    pair<handle_t, handle_t> key = make_pair(tip, tip);
                    
                    // If we already found the distance to itself, get it.
                    // Otherwise we'll do a Dijkstra from each tip (but only the
                    // ones we didn't do when looking at distinct tips).
                    auto tip_distance = get_or_compute_distance(key);
                    
                    if (tip_distance != numeric_limits<size_t>::max() &&
                        (tip_distance > furthest_distance || furthest_distance == numeric_limits<size_t>::max())) {
                        // If it's a new longer distance (but not infinite), between distinct tips, keep it
                        furthest_distance = tip_distance;
                        furthest = key;
                    }
                }

            }
            
            if (furthest_distance != numeric_limits<size_t>::max()) {
                // We found something!
                
#ifdef debug
                cerr << "Furthest Dijkstra shortest path between tips: " << furthest_distance << " bp" << endl;
#endif
                
                // Use the furthest-apart pair of tips as our telomeres
                add_telomeres(furthest.first, furthest.second);
                continue;
            } else {
#ifdef debug
                cerr << "No Dijkstra path found between any two tips." << endl;
#endif
            }
        }
        
        // Try even with disconnected tips before breaking into cycles. TODO:
        // Should we do this in preference to cycle breaking? We may eventually
        // want to try it both ways and compare.
        {
            if (component_tips[i].size() >= 2) {
                // TODO: For now we just pick two arbitrary tips in each component.
                vector<handle_t> tips{component_tips[i].begin(), component_tips[i].end()};
                add_telomeres(tips.front(), tips.back());
                continue;
            }
        }

        
        // Otherwise, we have no pair of tips. We have to
        // be cyclic, so find the biggest cycle (strongly connected component)
        // in the weakly connected component, pick an arbitrary node, and pick
        // the outward-facing ends of the node.
        
        {
            // What strongly connected components do we have?
            auto& strong_components = component_strong_components[i];
            
            assert(!strong_components.empty());
            
            // Find the largest in node size
            // Start out with no component selected
            size_t largest_component = 0;
            size_t largest_component_nodes = 0;
            for (size_t j = 0; j < strong_components.size(); j++) {
                // For each other strong component
                if (strong_components[j].size() > largest_component_nodes) {
                    // If it has more nodes, take it.
                    
                    if (strong_components[j].size() == 1) {
                        // Special check: Is this a real cycle? Nodes not in any
                        // cycle also get a strongly connected component size of
                        // 1.
                        
                        // Get the node we care about
                        handle_t member = graph.get_handle(*strong_components[j].begin(), false);
                        
                        // See if it cycles around.
                        // TODO: can we poll for the edge more cheaply?
                        bool saw_self = false;
                        graph.follow_edges(member, false, [&](const handle_t& next) {
                            if (next == member) {
                                saw_self = true;
                                return false;
                            }
                            return true;
                        });
                        
                        if (!saw_self) {
                            // This isn't really a 1-node cycle
#ifdef debug
                            cerr << "Component " << j << " is really a non-cyclic node and not a 1-node cycle" << endl;
#endif
                            continue;
                        }
                    }
                    
                    largest_component = j;
                    largest_component_nodes = strong_components[j].size();
                }
            }
            
            // Pick some arbitrary node in the largest component.
            // Also, assert we found one.
            assert(largest_component_nodes != 0);
            id_t break_node = *strong_components[largest_component].begin();
            
#ifdef debug
            cerr << "Elect to break at node " << break_node
                << " in largest strongly connected component with size " << largest_component_nodes << endl;
#endif
            
            // Break on it
            // Make sure to feed in telomeres facing out
            add_telomeres(graph.get_handle(break_node, true), graph.get_handle(break_node, false));
            continue;
        }
        
    }
    
    return make_pair(cactus_graph, telomeres);
}

VG cactus_to_vg(stCactusGraph* cactus_graph) {
    VG vg_graph;
    unordered_map<stCactusNode*, Node*> node_map;

    // keep track of mapping between nodes in graph
    function<Node*(stCactusNode*)> map_node = [&](stCactusNode* cac_node) -> Node* {
#ifdef debug
        cerr << "Node " << cac_node << " has object " << stCactusNode_getObject(cac_node)
            << " = " << *(id_t*)stCactusNode_getObject(cac_node) << endl;
#endif
        if (node_map.find(cac_node) == node_map.end()) {
            // Make sure to push IDs up by 1 because component numbering starts at 0    
            id_t cac_id = *(id_t*)stCactusNode_getObject(cac_node) + 1;
            Node* vg_node = vg_graph.create_node("N", cac_id + 1);
            node_map[cac_node] = vg_node;
            return vg_node;
        } else {
            return node_map[cac_node];
        }
    };

    // iterate nodes in cactus graph
    stCactusGraphNodeIt* node_it = stCactusGraphNodeIterator_construct(cactus_graph);
    stCactusNode* cac_node = stCactusGraphNodeIterator_getNext(node_it);
    
    for (; cac_node != NULL; cac_node = stCactusGraphNodeIterator_getNext(node_it)) {
        // iterate edges of node
        stCactusNodeEdgeEndIt edge_it = stCactusNode_getEdgeEndIt(cac_node);
        stCactusEdgeEnd* edge_end = stCactusNodeEdgeEndIt_getNext(&edge_it);
        for (; edge_end != NULL; edge_end = stCactusNodeEdgeEndIt_getNext(&edge_it)) {
            pair<stCactusNode*, stCactusNode*> cac_edge = make_pair(
                stCactusEdgeEnd_getNode(edge_end),
                stCactusEdgeEnd_getOtherNode(edge_end));
            // make a new vg edge
            Edge* vg_edge = vg_graph.create_edge(map_node(cac_edge.first),
                                                 map_node(cac_edge.second));
        }
    }
    stCactusGraphNodeIterator_destruct(node_it);

    return vg_graph;
}

VG cactusify(VG& graph) {
    if (graph.size() == 0) {
        return VG();
    }
    auto parts = handle_graph_to_cactus(graph, unordered_set<string>());
    stList_destruct(parts.second);
    stCactusGraph* cactus_graph = parts.first;
    VG out_graph = cactus_to_vg(cactus_graph);
    stCactusGraph_destruct(cactus_graph);
    return out_graph;
}

}
