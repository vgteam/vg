#include <unordered_set>
#include "bubbles.hpp"
#include "vg.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/weakly_connected_components.hpp"

extern "C" {
#include "sonLib.h"
#include "stCactusGraphs.h"
}

namespace vg {

using namespace std;
using namespace supbub;

SB_Input vg_to_sb_input(VG& graph){
    //cout << this->edge_count() << endl;
    SB_Input sbi;
    sbi.num_vertices = graph.edge_count();
    function<void(Edge*)> lambda = [&sbi](Edge* e){
        //cout << e->from() << " " << e->to() << endl;
        pair<id_t, id_t> dat = make_pair(e->from(), e->to() );
        sbi.edges.push_back(dat);
    };
    graph.for_each_edge(lambda);
    return sbi;
}

vector<pair<id_t, id_t> > get_superbubbles(SB_Input sbi){
    vector<pair<id_t, id_t> > ret;
    supbub::Graph sbg (sbi.num_vertices);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST superBubblesList{};
    supbub::DetectSuperBubble dsb;
    dsb.find(sbg, superBubblesList);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST::iterator it;
    for (it = superBubblesList.begin(); it != superBubblesList.end(); ++it) {
        ret.push_back(make_pair((*it).entrance, (*it).exit));
    }
    return ret;
}

vector<pair<id_t, id_t> > get_superbubbles(VG& graph){
    vector<pair<id_t, id_t> > ret;
    supbub::Graph sbg (graph.max_node_id() + 1);
    //load up the sbgraph with edges
    function<void(Edge*)> lambda = [&sbg](Edge* e){
#ifdef debug
        cout << e->from() << " " << e->to() << endl;
#endif
        sbg.addEdge(e->from(), e->to());
    };

    graph.for_each_edge(lambda);

    supbub::DetectSuperBubble::SUPERBUBBLE_LIST superBubblesList{};

    supbub::DetectSuperBubble dsb;
    dsb.find(sbg, superBubblesList);
    supbub::DetectSuperBubble::SUPERBUBBLE_LIST::iterator it;
    for (it = superBubblesList.begin(); it != superBubblesList.end(); ++it) {
        ret.push_back(make_pair((*it).entrance, (*it).exit));
    }
    return ret;
}
// check for conflict (duplicate nodes and edges) occurs within add_* functions

map<pair<id_t, id_t>, vector<id_t> > superbubbles(VG& graph) {
    map<pair<id_t, id_t>, vector<id_t> > bubbles;
    // flip doubly reversed edges
    graph.flip_doubly_reversed_edges();
    // ensure we're sorted
    algorithms::sort(&graph);
    // if we have a DAG, then we can find all the nodes in each superbubble
    // in constant time as they lie in the range between the entry and exit node
    auto supbubs = get_superbubbles(graph);
    //     hash_map<Node*, int> node_index;
    for (auto& bub : supbubs) {
        auto start = graph.node_index[graph.get_node(bub.first)];
        auto end = graph.node_index[graph.get_node(bub.second)];
        // get the nodes in the range
        auto& b = bubbles[bub];
        for (int i = start; i <= end; ++i) {
            auto node = graph.graph.mutable_node(i);
            if (i == start || i == end
                || !graph.is_head_node(node)
                && !graph.is_tail_node(node)) {
                b.push_back(node->id());
            }
        }
    }
    return bubbles;
}

// Cactus


// Step 1) Get undirected adjacency connected components of VG *sides*
// (note we can't put NodeSides in header, so leave this function static)
typedef set<NodeSide> SideSet;
typedef map<NodeSide, int> Side2Component;
static void compute_side_components(VG& graph, 
                                    vector<SideSet>& components,
                                    Side2Component& side_to_component) {

    // update components datastructures with a side.  if comp is -1
    // create new component
    function<void(NodeSide, int)> update_component = [&] (NodeSide side, int comp) {
        if (comp == -1) {
            components.push_back(SideSet());
            comp = components.size() - 1;
        }
        SideSet& component = components[comp];
        assert(side_to_component.find(side) == side_to_component.end());
        side_to_component[side] = comp;
        component.insert(side);
    };
    
    // create a connected component of a node side and everything adjacent
    // (if not added already)
    function<void(NodeSide)> add_node_side = [&](NodeSide side) {
        queue<NodeSide> q;
        q.push(side);
        while (!q.empty()) {
            NodeSide cur_side = q.front();
            q.pop();
            if (side_to_component.find(cur_side) == side_to_component.end()) {
                update_component(cur_side, cur_side == side ? -1 : components.size() - 1);
                // visit all adjacent sides
                set<NodeSide> adj_sides = graph.sides_of(cur_side);
                for (auto adj_side : adj_sides) {
                    q.push(adj_side);
                }
            }
        }
    };

    graph.for_each_node([&](Node* n) {
            add_node_side(NodeSide(n->id(), false));
            add_node_side(NodeSide(n->id(), true));
        });
}

struct CactusSide {
    int64_t node;
    bool is_end;
};

void* mergeNodeObjects(void* a, void* b) {
    // One of the objects is going to get returned, and the other one is going
    // to get freed (since the graph is supposed to won them all).
    id_t* to_return;
    id_t* to_free;
    
    if(*(id_t*)a < *(id_t*)b) {
        to_return = (id_t*)a;
        to_free = (id_t*)b;
    } else {
        to_free = (id_t*)a;
        to_return = (id_t*)b;
    }
    
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

#define debug
// Step 2) Make a Cactus Graph. Returns the graph and a list of paired
// cactusEdgeEnd telomeres, one after the other. Both members of the return
// value must be destroyed.
pair<stCactusGraph*, stList*> vg_to_cactus(VG& graph) {

    // in a cactus graph, every node is an adjacency component.
    // every edge is a *vg* node connecting the component

    // start by identifying adjacency components
    vector<SideSet> components;
    Side2Component side_to_component;
    compute_side_components(graph, components, side_to_component);

    // map cactus nodes back to components
    vector<stCactusNode*> cactus_nodes(components.size());
    
    // create cactus graph
    stCactusGraph* cactus_graph = stCactusGraph_construct2(free, free);
    
    // copy each component into cactus node
    for (int i = 0; i < components.size(); ++i) {
#ifdef debug
        cout << "Creating cactus node for component " << i << " with size " << components[i].size() << endl;
        for (auto& content : components[i]) {
            cout << "\tContains: " << content << endl;
        }
#endif
        id_t* cactus_node_id = (id_t*)malloc(sizeof(id_t));
        *cactus_node_id = i;
        cactus_nodes[i] = stCactusNode_construct(cactus_graph, cactus_node_id);
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
            NodeSide other_side(side.node, !side.is_end);
            // And what component it is in
            int j = side_to_component[other_side];
            
            if (!edge_ends.count(side.node)) {
                // If we haven't made the Cactus edge for this graph node yet
            
                // afraid to try to get C++ NodeSide class into C, so we copy to
                // equivalent struct
                CactusSide* cac_side1 = (CactusSide*)malloc(sizeof(CactusSide));
                CactusSide* cac_side2 = (CactusSide*)malloc(sizeof(CactusSide));
                cac_side1->node = side.node;
                cac_side1->is_end = side.is_end;
                cac_side2->node = other_side.node;
                cac_side2->is_end = other_side.is_end;
#ifdef debug
                cout << "Creating cactus edge for sides " << side << " -- " << other_side << ": "
                     << i << " -> " << j << endl;
#endif

                // We get the cactusEdgeEnd corresponding to the side stored in side.
                // This may be either the left (if that NodeSide is a start), or the right (if that NodeSide is an end).
                stCactusEdgeEnd* cactus_edge = stCactusEdgeEnd_construct(
                    cactus_graph,
                    cactus_nodes[i],
                    cactus_nodes[j],
                    cac_side1,
                    cac_side2);
                // Save the cactusEdgeEnd for the left side of the node.
                edge_ends[side.node] = side.is_end ? stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge) : cactus_edge;
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
        cerr << "Found head " << graph.get_id(head) << " in component " << node_to_component[graph.get_id(head)] << endl;
    }
    for (auto& tail : all_tails) {
        component_tips[node_to_component[graph.get_id(tail)]].insert(graph.flip(tail));
        cerr << "Found tail " << graph.get_id(tail) << " in component " << node_to_component[graph.get_id(tail)] << endl;
    }
    
    // Assign path names to components
    vector<vector<string>> component_paths(weak_components.size());
    // Also get the path length.
    unordered_map<string, size_t> path_length;
    
    graph.paths.for_each_name([&](const std::string& name) {
        // For every path
        auto& path_mappings = graph.paths.get_path(name);
        
        // Save the path under the component
        auto component = node_to_component[path_mappings.front().position().node_id()];
        component_paths[component].push_back(name);
        
        cerr << "Path " << name << " belongs to component " << component << endl;
        
        for (auto& mapping : graph.paths.get_path(name)) {
            // Total up the length. We could use from length on the mapping, but
            // sometimes (like in the tests) the mapping edits haven't been
            // populated.
            path_length[name] += graph.get_length(graph.get_handle(mapping.position().node_id(), false));
        }
        
        cerr << "\tPath " << name << " has length " << path_length[name] << endl;
    });
    
    
    // OK, now we need to fill in the telomeres list with two telomeres per
    // component.
    
    // This function adds as telomeres a pair of inward-facing handles.
    auto add_telomeres = [&](const handle_t& left, const handle_t& right) {
        cerr << "Selected " << graph.get_id(left) << " " << graph.get_is_reverse(left) << " and "
            << graph.get_id(right) << " " << graph.get_is_reverse(right) << " as tips" << endl;
        
        stCactusEdgeEnd* end1 = edge_ends[graph.get_id(left)];
        // We need to add the interior side of the node, and our handle is reading inwards.
        // If we're reverse, add the node's local left. Otherwise, add the node's local right.
        stList_append(telomeres, graph.get_is_reverse(left) ? end1 : stCactusEdgeEnd_getOtherEdgeEnd(end1));
        stCactusEdgeEnd* end2 = edge_ends[graph.get_id(right)];
        stList_append(telomeres, graph.get_is_reverse(right) ? end2 : stCactusEdgeEnd_getOtherEdgeEnd(end2));
    };
    
    for (size_t i = 0; i < weak_components.size(); i++) {
        // For each weakly connected component
        
        // First priority is two tips at the ends of the longest path that has
        // two tips at its ends.
        {
            string longest_path;
            // This is going to hold inward-facing handles to the tips we find.
            pair<handle_t, handle_t> longest_path_tips;
            size_t longest_path_length = 0;
            cerr << "Consider " << component_paths[i].size() << " paths for component " << i << endl;
            for (auto& path_name : component_paths[i]) {
                // Look at each path
                auto& path_mappings = graph.paths.get_path(path_name);
                
                cerr << "\tPath " << path_name << " has " << path_mappings.size() << " mappings" << endl;
                
                // See if I can get two tips on its ends.
                // Get the inward-facing start and end handles.
                handle_t path_start = graph.get_handle(path_mappings.front().position().node_id(),
                    path_mappings.front().position().is_reverse());
                handle_t path_end = graph.get_handle(path_mappings.back().position().node_id(),
                    !path_mappings.back().position().is_reverse());
                
                if (component_tips[i].count(path_start) && component_tips[i].count(path_end)) {
                    // This path ends in two tips so we can consider it
                    
                    if (path_length[path_name] > longest_path_length) {
                        // This is our new longest path between tips.
                        longest_path = path_name;
                        longest_path_tips = make_pair(path_start, path_end);
                        longest_path_length = path_length[path_name];
                        cerr << "\t\tNew longest path!" << endl;
                    } else {
                        cerr << "\t\tPath length of " << path_length[path_name] << " not longer than path "
                            << longest_path << " with " << longest_path_length << endl;
                    } 
                    
                } else {
                    cerr << "\t\tPath " << path_name << " does not start and end with tips" << endl;
                }
            }
            
            if (!longest_path.empty()) {
                cerr << "Longest tip path is " << longest_path << endl;
                // We found something!
                add_telomeres(longest_path_tips.first, longest_path_tips.second);
                // Work on the next component.
                continue;
            } else {
                cerr << "No longest path found between tips" << endl;
            }
        }
        
        // Second priority is Dijkstra traversal tip selection
        
        // TODO: Otherwise, compute tip reachability and distance by Dijkstra's
        // Algorithm and pick the pair of reachable tips with the longest shortest
        // paths
        
        // TODO: Otherwise, we have to be cyclic, so find the biggest cycle
        // (strongly connected component) in the weakly connected component, pick an
        // arbitrary node, and pick the outward-facing ends of the node.

        // TODO: For now we just pick two arbitrary tips in each component.
        vector<handle_t> tips{component_tips[i].begin(), component_tips[i].end()};
        add_telomeres(tips.front(), tips.back());
        
    }
    
    return make_pair(cactus_graph, telomeres);
}
#undef debug

// fill in the "acyclic" and "contents" field of a bubble by doing a depth first search
// between its bounding sides (start and end)
static void fill_ultrabubble_contents(VG& graph, Bubble& bubble) {

    // orient out from source
    NodeTraversal source(graph.get_node(bubble.start.node), !bubble.start.is_end);
    vector<NodeTraversal> sources = {source};
    // but never walk "through" sink in any direction
    // note we treat the source in the opposite direction as a sink here to prevent leaving
    // the bubble in a loop
    id_t sink_id = bubble.end.node;
    unordered_set<NodeTraversal> sinks = {NodeTraversal(graph.get_node(bubble.start.node), bubble.start.is_end),
                                          NodeTraversal(graph.get_node(sink_id), !bubble.start.is_end),
                                          NodeTraversal(graph.get_node(sink_id), bubble.start.is_end)};

    // remember unique node ids we've visited 
    unordered_set<id_t> contents_set;
    
    // the acyclic logic derived from vg::is_acyclic()
    // but changed to make sure we only ever touch the ends of the
    // source and sink (never loop over body of the node)
    // and to look for directed (as opposed to bidirected) cycles
    unordered_set<vg::id_t> seen;
    bubble.dag = true;
    
    graph.dfs([&](NodeTraversal trav) {

            // self loops on bubble end points considered degenerate
            if (trav.node->id() != source.node->id() && trav.node->id() != sink_id &&
                graph.is_self_looping(trav.node)) {
                bubble.dag = false;
            }
            // don't step past the sink once we've reached it
            if (sinks.count(trav) == false) {
                for (auto& next : graph.travs_from(trav)) {
                    // filter out self loop on bubble endpoints like above
                    if (trav.node->id() != source.node->id() && next.node->id() != trav.node->id() &&
                        seen.count(next.node->id())) {
                        bubble.dag = false;
                        break;
                    }
                }
            }
            if (bubble.dag) {
                seen.insert(trav.node->id());
            }

            contents_set.insert(trav.node->id());
        },
        [&](NodeTraversal trav) {
            seen.erase(trav.node->id());
        },
        &sources,
        &sinks);

    bubble.contents.clear();
    bubble.contents.insert(bubble.contents.begin(), contents_set.begin(), contents_set.end());
}

// cactus C data to C++ tree interface (ultrabubble as added to child_list of out_node)
// filling in the internace bubble nodes as well as acyclicity using the original graph
static void ultrabubble_recurse(VG& graph, stList* chains_list,
                                NodeSide side1, NodeSide side2, BubbleTree::Node* out_node) {

    // add the Tree node
    out_node->v.start = side1;
    out_node->v.end = side2;

    // only time this won't be true, is for the dummy root node
    if (side1.node != 0) {
        assert(side2.node != 0);
        fill_ultrabubble_contents(graph, out_node->v);

    } 
        
    int chain_offset = 0;
    // for each nested chain
    for (int64_t i = 0; i < stList_length(chains_list); i++) {

        stList* cactus_chain = (stList*)stList_get(chains_list, i);
            
        // for each ultra bubble in the chain
        for (int64_t j = 0; j < stList_length(cactus_chain); j++) {

            // add a new chain offset to our output list, remembering where chain started
            if (j == 0) {
                out_node->v.chain_offsets.push_back(chain_offset);                    
            }
            // add a node to our output tree
            BubbleTree::Node* new_node = new BubbleTree::Node();
            out_node->children.push_back(new_node);
            new_node->parent = out_node;
            
            stSnarl* child_bubble = (stSnarl*)stList_get(cactus_chain, j);

            // scrape the vg coordinate information out of the cactus ends where we stuck
            // it during cactus construction
            CactusSide* cac_child_side1 = (CactusSide*)stCactusEdgeEnd_getObject(child_bubble->edgeEnd1);
            CactusSide* cac_child_side2 = (CactusSide*)stCactusEdgeEnd_getObject(child_bubble->edgeEnd2);
            NodeSide child_side1(cac_child_side1->node, cac_child_side1->is_end);
            NodeSide child_side2(cac_child_side2->node, cac_child_side2->is_end);

            // try to keep bubble sides sorted to be more consistent with superbubbles
            if (child_side2 < child_side1) {
                std::swap(child_side1, child_side2);
            }
                
            ultrabubble_recurse(graph, child_bubble->chains, child_side1, child_side2, new_node);
                
            ++chain_offset;
        }
    }
}

BubbleTree* ultrabubble_tree(VG& graph) {
    if (graph.size() == 0) {
        return new BubbleTree(new BubbleTree::Node());
    }
    // convert to cactus
    pair<stCactusGraph*, stList*> cac_pair = vg_to_cactus(graph);
    stCactusGraph* cactus_graph = cac_pair.first;
    stList* telomeres = cac_pair.second;

    BubbleTree* out_tree = new BubbleTree(new BubbleTree::Node());

    // get the snarl decomposition as a C struct
    stSnarlDecomposition *snarls = stCactusGraph_getSnarlDecomposition(cactus_graph, telomeres);
    
    // Get a non-owning pointer to the list of chains (which are themselves lists of snarls).
    stList* cactus_chains_list = snarls->topLevelChains;
    
    // copy back to our C++ tree interface (ultrabubble as added to child_list of out_node)
    // in to new ultrabubble code, we no longer have tree root.  instead, we shoehorn
    // dummy root onto tree just to preserve old interface for now.
    ultrabubble_recurse(graph, cactus_chains_list, NodeSide(), NodeSide(), out_tree->root);
    
    // Free the decomposition
    stSnarlDecomposition_destruct(snarls);
    
    // Free the telomeres
    stList_destruct(telomeres);

    // free the cactus graph
    stCactusGraph_destruct(cactus_graph);

    return out_tree;
}

map<pair<id_t, id_t>, vector<id_t> > ultrabubbles(VG& graph) {

    algorithms::sort(&graph);
    map<pair<id_t, id_t>, vector<id_t> > output;

    BubbleTree* bubble_tree = ultrabubble_tree(graph);

    bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {
            // cut root to be consistent with superbubbles()
            if (node != bubble_tree->root) {
                Bubble& bubble = node->v;
                // sort nodes to be consistent with superbubbles
                sort(bubble.contents.begin(), bubble.contents.end());
                output[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
            }
        });

    delete bubble_tree;
    return output;
}

VG cactus_to_vg(stCactusGraph* cactus_graph) {
    VG vg_graph;
    unordered_map<stCactusNode*, Node*> node_map;

    // keep track of mapping between nodes in graph
    function<Node*(stCactusNode*)> map_node = [&](stCactusNode* cac_node) -> Node* {
        if (node_map.find(cac_node) == node_map.end()) {
            id_t cac_id = *(id_t*)stCactusNode_getObject(cac_node);
            Node* vg_node = vg_graph.create_node("N", cac_id);
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
            // how to use?
            CactusSide* cac_side = (CactusSide*)stCactusEdgeEnd_getObject(
                stCactusEdgeEnd_getLink(edge_end));
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
    auto parts = vg_to_cactus(graph);
    stList_destruct(parts.second);
    stCactusGraph* cactus_graph = parts.first;
    VG out_graph = cactus_to_vg(cactus_graph);
    stCactusGraph_destruct(cactus_graph);
    return out_graph;
}

}
