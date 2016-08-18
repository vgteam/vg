#include <unordered_set>
#include "bubbles.hpp"
#include "vg.hpp"

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
    graph.sort();
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
static void compute_side_components(VG& graph, pair<NodeSide, NodeSide> source_sink,
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

    // make component for source--sink edge
    NodeSide first_side(source_sink.first.node, !source_sink.first.is_end);
    NodeSide last_side(source_sink.second.node, !source_sink.second.is_end);
    update_component(first_side, -1);
    update_component(last_side, 0);
    
    // create a connected component of a node side and everything adjancet
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

// Step 2) Make a Cactus Graph
pair<stCactusGraph*, stCactusNode*> vg_to_cactus(VG& graph, pair<NodeSide, NodeSide> source_sink) {

    // in a cactus graph, every node is an adjacency component.
    // every edge is a *vg* node connecting the component

    // start by identifying adjacency components
    vector<SideSet> components;
    Side2Component side_to_component;
    compute_side_components(graph, source_sink, components, side_to_component);

    // map cactus nodes back to components
    vector<stCactusNode*> cactus_nodes(components.size());

    // create cactus graph
    stCactusGraph* cactus_graph = stCactusGraph_construct2(free, free);
    
    // copy each component into cactus node
    for (int i = 0; i < components.size(); ++i) {
#ifdef debug
        cout << "Creating cactus node for compontnet " << i << " with size " << components[i].size() << endl;
#endif
        id_t* cactus_node_id = (id_t*)malloc(sizeof(id_t));
        *cactus_node_id = i;
        cactus_nodes[i] = stCactusNode_construct(cactus_graph, cactus_node_id);
    }

    // make edge for each vg node connecting two components
    // they are undirected, so we use this set to keep track
    unordered_set<id_t> created_edges;
    for (int i = 0; i < components.size(); ++i) {
        for (auto side : components[i]) {
            NodeSide other_side(side.node, !side.is_end);
            int j = side_to_component[other_side];
            if (created_edges.find(side.node) == created_edges.end()) {
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
                stCactusEdgeEnd* cactus_edge = stCactusEdgeEnd_construct(
                    cactus_graph,
                    cactus_nodes[i],
                    cactus_nodes[j],
                    cac_side1,
                    cac_side2);
                created_edges.insert(side.node);
            }
        }
    }

    // collapse 3-edge connected components to make actual cactus graph. 
    stCactusGraph_collapseToCactus(
        cactus_graph, mergeNodeObjects, cactus_nodes[0]);

    // don't worry about bridges
    stCactusGraph_collapseBridges(
        cactus_graph, cactus_nodes[0], mergeNodeObjects);

    return make_pair(cactus_graph, cactus_nodes[0]);

}

pair<NodeSide, NodeSide> get_cactus_source_sink(VG& graph) {

    auto source_ids = graph.head_nodes();
    auto sink_ids = graph.tail_nodes();

    // todo : relax?
    assert(!source_ids.empty() && !sink_ids.empty());

    // important : graph must be sorted!

    // using these ranks requires graph.sort() to have been run
    // to make sense
    int source_rank = graph.size();
    for (auto n : source_ids) {
      source_rank = min(source_rank, graph.node_rank(n));
    }
    int sink_rank = 0;
    for (auto n : sink_ids) {
      sink_rank = max(sink_rank, graph.node_rank(n));
    }

    // go back from rank to node
    id_t source_id = graph.graph.node(source_rank).id();
    id_t sink_id = graph.graph.node(sink_rank).id();

    assert(source_id != sink_id);

    return make_pair(NodeSide(source_id, true), NodeSide(sink_id, false));
}

pair<NodeSide, NodeSide> get_cactus_source_sink(VG& graph, const string& path_name, int steps) {

    list<Mapping>& path = graph.paths.get_path(path_name);
    
    NodeSide source(path.front().position().node_id(), true);
    NodeSide sink(path.back().position().node_id(), false);

    // search context of path ends for actual head and tail nodes
    function<pair<id_t, id_t>(id_t)> find_endpoints = [&](id_t node_id) -> pair<id_t, id_t> {
        pair<id_t, id_t> found_ends(graph.is_head_node(node_id) ? node_id : 0,
                                    graph.is_tail_node(node_id) ? node_id : 0);
        VG context;
        graph.nonoverlapping_node_context_without_paths(graph.get_node(node_id), context);
        graph.expand_context(context, steps, false);
        context.for_each_node([&](Node* node) {
                if (found_ends.first == 0 && graph.is_head_node(node->id())) {
                    // todo: pick best from many?
                    found_ends.first = node->id();
                }
                if (found_ends.second == 0 && graph.is_tail_node(node->id())) {
                    found_ends.second = node->id();
                }
            });
        return found_ends;
    };

    auto source_ends = find_endpoints(source.node);
    auto sink_ends = find_endpoints(sink.node);

    // favour  head/tail or tail/head matchup
    if (source_ends.first > 0 && sink_ends.second > 0) {
        source = NodeSide(source_ends.first, true);
        sink = NodeSide(sink_ends.second, false);
    } else if (source_ends.second > 0 && sink_ends.first > 0) {
        source = NodeSide(source_ends.second, false);
        sink = NodeSide(sink_ends.first, true);
    } // settle for head / head or tail / tail (to what results?)
    else if (source_ends.first > 0 && sink_ends.first > 0) {
        source = NodeSide(source_ends.first, true);
        sink = NodeSide(sink_ends.first, true);
    } else if (source_ends.second > 0 && sink_ends.second > 0) {
        source = NodeSide(source_ends.second, false);
        sink = NodeSide(sink_ends.second, false);
    }

    assert((source.is_end && graph.is_head_node(source.node)) ||
           (!source.is_end && graph.is_tail_node(source.node)));
    assert((sink.is_end && graph.is_head_node(sink.node)) ||
           (!sink.is_end && graph.is_tail_node(sink.node)));
    
    return make_pair(source, sink);
}

void remove_heads_and_tails(VG& graph, pair<NodeSide, NodeSide> source_sink) {
    // collect all heads and tails
    vector<id_t> heads_tails;
    graph.for_each_node([&](Node* node) {
            if (graph.is_head_node(node->id()) || graph.is_tail_node(node->id())) {
                heads_tails.push_back(node->id());
            }
        });

    // iteratively prune until the source and sink are the only heads and tails.
    vector<id_t> next;
    while (!heads_tails.empty()) {
        for (auto node_id : heads_tails) {
            if (graph.has_node(node_id) && // node can linked to from multiple source/sinks
                (graph.is_head_node(node_id) || graph.is_tail_node(node_id)) &&
                node_id != source_sink.first.node && node_id != source_sink.second.node) {
                for (auto adj : graph.sides_context(node_id)) {                
                    next.push_back(adj.first.node);
                }
                graph.destroy_node(node_id);
            }
        }
        swap(heads_tails, next);
        next.clear();
    }
}

BubbleTree cactusbubble_tree(VG& graph, pair<NodeSide, NodeSide> source_sink) {
    
    // convert to cactus
    pair<stCactusGraph*, stCactusNode*> cac_pair = vg_to_cactus(graph, source_sink);
    stCactusGraph* cactus_graph = cac_pair.first;
    stCactusNode* root_node = cac_pair.second;

    BubbleTree out_tree(new BubbleTree::Node());
        
    // recursive function to walk through the cycles in the cactus graph
    function<vector<id_t>(stCactusNode*, BubbleTree::Node*)> cactus_recurse = [&](stCactusNode* cactus_node,
                                                                                  BubbleTree::Node* out_node) {

        // walk through the edges incident with the node and find those that are at the
        // start of a chain
        stCactusNodeEdgeEndIt edge_it = stCactusNode_getEdgeEndIt(cactus_node);
        stCactusEdgeEnd* cactus_edge_end;
        // A list of all the edges in the superbubble.
        vector<id_t> edgesVisited;
        
        while ((cactus_edge_end = stCactusNodeEdgeEndIt_getNext(&edge_it)) != NULL) {
            
            if (stCactusEdgeEnd_isChainEnd(cactus_edge_end) &&
                stCactusEdgeEnd_getLinkOrientation(cactus_edge_end)) {
                // Add the edge to those visited
                CactusSide* side = (CactusSide*)stCactusEdgeEnd_getObject(cactus_edge_end);
                edgesVisited.push_back(side->node);
                
                // Walk the chain starting from the first edge in the chain
                assert(stCactusEdgeEnd_getNode(cactus_edge_end) == cactus_node);
                stCactusEdgeEnd* cactus_edge_end2 = stCactusEdgeEnd_getLink(
                    stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge_end));
                
                while(cactus_edge_end != cactus_edge_end2) {

                    // Report super bubble!!!!!!!
                    CactusSide* side1 = (CactusSide*)stCactusEdgeEnd_getObject(cactus_edge_end2);
                    CactusSide* side2 = (CactusSide*)stCactusEdgeEnd_getObject(
                        stCactusEdgeEnd_getLink(cactus_edge_end2));

                    // Add the next edge in the chain to the edges visited
                    edgesVisited.push_back(side2->node);
                    
                    assert(stCactusEdgeEnd_getNode(cactus_edge_end2) != cactus_node);

                    BubbleTree::Node* new_node = new BubbleTree::Node();
                    out_node->children.push_back(new_node);

                    // This recursive explores the cycles in the subgraph rooted at next_cactus_node.
                    new_node->v.start = NodeSide(side1->node, side1->is_end);
                    new_node->v.end = NodeSide(side2->node, side2->is_end);
                    new_node->v.contents = cactus_recurse(stCactusEdgeEnd_getNode(cactus_edge_end2), new_node);
                    new_node->v.contents.insert(new_node->v.contents.begin(), side1->node);
                    new_node->v.contents.push_back(side2->node);

                    // Add all the edges in the super bubble to the total list of edges
                    edgesVisited.insert(edgesVisited.begin(),
                                        out_node->v.contents.begin(),
                                        out_node->v.contents.end());

                    // Get the next link
                    cactus_edge_end2 = stCactusEdgeEnd_getLink(
                        stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge_end2));
                }
            }
        }
        return edgesVisited;
    };
    
    cactus_recurse(root_node, out_tree.root);

    stCactusGraph_destruct(cactus_graph);

    // compute root as special case (todo: smarten up)
    //  these were used for source/sink in compute_side_components() above
    out_tree.root->v.start = source_sink.first;
    out_tree.root->v.end = source_sink.second;
    
    return out_tree;
}

void bubble_up_bubbles(BubbleTree& bubble_tree) {
    bubble_tree.for_each_postorder([&](BubbleTree::Node* node) {
            for (auto& c: node->children) {
                node->v.contents.insert(node->v.contents.end(),
                                       c->v.contents.begin(), c->v.contents.end());
                
            }
            if (!node->children.empty()) {
                // todo : uniqueness shouldnt be issue -- debug while sober
                set<id_t> unique_set{node->v.contents.begin(), node->v.contents.end()};
                node->v.contents.clear();
                node->v.contents.insert(node->v.contents.begin(), unique_set.begin(), unique_set.end());
            }
        });
}

void bubble_down_bubbles(BubbleTree& bubble_tree) {
    bubble_tree.for_each_preorder([&](BubbleTree::Node* node) {
            set<id_t> unique_set{node->v.contents.begin(), node->v.contents.end()};
            if (node == bubble_tree.root) {
                node->v.contents = vector<id_t>{*unique_set.begin(), *unique_set.end()};
            } else {
                set<id_t> child_set;
                for (auto& c: node->children) {
                    child_set.insert(c->v.contents.begin(), c->v.contents.end());                
                }
                auto it = set_difference(unique_set.begin(), unique_set.end(), child_set.begin(), child_set.end(),
                               node->v.contents.begin());
                node->v.contents.resize(it - node->v.contents.begin());
            }
        });
}

map<pair<id_t, id_t>, vector<id_t> > cactusbubbles(VG& graph) {

    graph.sort();
    map<pair<id_t, id_t>, vector<id_t> > output;

    // get endpoints
    pair<NodeSide, NodeSide> source_sink = get_cactus_source_sink(graph);

    BubbleTree bubble_tree = cactusbubble_tree(graph, source_sink);
    bubble_up_bubbles(bubble_tree);

    bubble_tree.for_each_preorder([&](BubbleTree::Node* node) {
            Bubble& bubble = node->v;
            // cut root to be consistent with superbubbles()
            if (bubble.start != bubble_tree.root->v.start ||
                bubble.end != bubble_tree.root->v.end) {               
                // sort nodes to be consistent with superbubbles
                sort(bubble.contents.begin(), bubble.contents.end());
                output[make_pair(bubble.start.node, bubble.end.node)] = bubble.contents;
            }
        });

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
    pair<NodeSide, NodeSide> source_sink = get_cactus_source_sink(graph);
    stCactusGraph* cactus_graph = vg_to_cactus(graph, source_sink).first;
    VG out_graph = cactus_to_vg(cactus_graph);
    stCactusGraph_destruct(cactus_graph);
    return out_graph;
}

}
