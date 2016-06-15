#include <unordered_set>
#include "bubbles.hpp"
#include "vg.hpp"

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
            b.push_back(graph.graph.node(i).id());
        }
    }
    return bubbles;
}

// Cactus
// Step 1) Get undirected adjacency connected components of VG *sides*
// (note we can't put NodeSides in header, so leave this function static)
typedef set<NodeSide> SideSet;
typedef map<NodeSide, int> Side2Component;
static void compute_side_components(VG& graph, vector<SideSet>& components,
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
    NodeSide first_side(graph.min_node_id(), false);
    NodeSide last_side(graph.max_node_id(), true);
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
    return *(id_t*)a < *(id_t*)b ? a : b; 
}

// Step 2) Make a Cactus Graph
pair<stCactusGraph*, stCactusNode*> vg_to_cactus(VG& graph) {

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

map<pair<id_t, id_t>, vector<id_t> > cactusbubbles(VG& graph) {

    // convert to cactus
    pair<stCactusGraph*, stCactusNode*> cac_pair = vg_to_cactus(graph);
    stCactusGraph* cactus_graph = cac_pair.first;
    stCactusNode* root_node = cac_pair.second;

    map<pair<id_t, id_t>, vector<id_t> > output;
        
    // recursive function to walk through the cycles in the cactus graph
    function<vector<id_t>(stCactusNode*)> cactusRecurse = [&](stCactusNode *cactus_node) {

        // walk through the edges incident with the node and find those that are at the
        // start of a chain
        stCactusNodeEdgeEndIt edge_it = stCactusNode_getEdgeEndIt(cactus_node);
        stCactusEdgeEnd *cactus_edge_end;
        vector<id_t> edgesVisited; // A list of all the edges in the superbubble.
        while(cactus_edge_end = stCactusNodeEdgeEndIt_getNext(&edge_it)) {
            if(stCactusEdgeEnd_isChainEnd(cactus_edge_end) && stCactusEdgeEnd_getLinkOrientation(cactus_edge_end))  {
                // Add the edge to those visited
                edgesVisited.push_back(((CactusSide*)stCactusEdgeEnd_getObject(cactus_edge_end))->node);
                
                //Walk the chain starting from the first edge in the chain
                assert(stCactusEdgeEnd_getNode(cactus_edge_end) == cactus_node);
                stCactusEdgeEnd *cactus_edge_end2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge_end));
                while(cactus_edge_end != cactus_edge_end2) {
                    // Add the next edge in the chain to the edges visited
                    edgesVisited.push_back(((CactusSide*)stCactusEdgeEnd_getObject(cactus_edge_end2))->node);
                    
                    assert(stCactusEdgeEnd_getNode(cactus_edge_end2) != cactus_node);

                    // Report super bubble!!!!!!!
                    CactusSide* side1 = (CactusSide*)stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(cactus_edge_end2));
                    CactusSide* side2 = (CactusSide*)stCactusEdgeEnd_getObject(cactus_edge_end2);

                    // This recursive explores the cycles in the subgraph rooted at next_cactus_node.
                    output[make_pair(side1->node, side2->node)] = cactusRecurse(stCactusEdgeEnd_getNode(cactus_edge_end2));
                    // Add all the edges in the super bubble to the total list of edges
                    edgesVisited.insert(edgesVisited.begin(),
                                        output[make_pair(side1->node, side2->node)].begin(),
                                        output[make_pair(side1->node, side2->node)].end());

                    // Get the next link
                    cactus_edge_end2 = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(cactus_edge_end2));
                }
            }
        }
        return edgesVisited;
    };
    
    cactusRecurse(root_node);

    return output;
}
    

}
