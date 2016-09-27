#ifndef BUBBLES_H

#define BUBBLES_H

#include <vector>
#include <map>

#include "types.hpp"
#include "utility.hpp"
#include "nodeside.hpp"

#include "DetectSuperBubble.hpp"
extern "C" {
    typedef struct _stCactusGraph stCactusGraph;
    typedef struct _stCactusNode stCactusNode;;
}

using namespace std;

namespace vg {

class VG;

// Consolidate bubble finding code here for both cactus and superbubbles
// to keep vg class size from getting even more out of hand

// SUPERBUBBLES

struct SB_Input{
    int num_vertices;
    vector<pair<id_t, id_t> > edges;
};

SB_Input vg_to_sb_input(VG& graph);
vector<pair<id_t, id_t> > get_superbubbles(SB_Input sbi);
vector<pair<id_t, id_t> > get_superbubbles(VG& graph);
map<pair<id_t, id_t>, vector<id_t> > superbubbles(VG& graph);

// CACTUS BUBBLES
// todo : refactor harmonize interface with stuff in deconstruct
// and superbubbles in general
struct Bubble {
    NodeSide start;
    NodeSide end;
    vector<id_t> contents;
    // cactus now gives us chaining information, stick here for now
    // so chain_offsets[i]-chain_offsets[i+1] mark the range
    // of children in chain i.  existing code that doesn't use
    // chains will be unaffected. 
    vector<int> chain_offsets;   
};

typedef Tree<Bubble> BubbleTree;

// Convert VG to Cactus Graph
// Notes:
//  - source and sink of vg graph obtained by min/max node id
//    so graph should be sorted
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns "root" node as well as graph
pair<stCactusGraph*, stCactusNode*> vg_to_cactus(VG& graph, pair<NodeSide, NodeSide> source_sink);

// Get source and sink nodes, relying on node ranks of sorted graph
pair<NodeSide, NodeSide> get_cactus_source_sink(VG& graph);

// Get source and sink nodes from path endpoints
pair<NodeSide, NodeSide> get_cactus_source_sink(VG& graph, const string& path_name, int steps = 5);

// Remove all head and tail nodes that aren't a source or sink
void remove_heads_and_tails(VG& graph, pair<NodeSide, NodeSide> source_sink);

// Return the hierchical cactus decomposition
// Input graph must be sorted!
BubbleTree cactusbubble_tree(VG& graph, pair<NodeSide, NodeSide> source_sink);

// Use the new ultrabubble interface in cactus. will replace above once tested
BubbleTree cactusbubble_tree_new(VG& graph, pair<NodeSide, NodeSide> source_sink);

// By default, bubble X's contents array doesn't have all the nodes
// in its children's contents (ie each node only stored in lowest bubble
// that contains it).  This function will store each node in every
// bubble that contains it (like superbubbles)
// Note: contents wont be in any particular order after
void bubble_up_bubbles(BubbleTree& bubble_tree);

// Revert the results of bubble_up_bubbles by subtracting
// child nodesets from parent's.  (also a temp workaround
// to fix up tree returned by cactusbubble_tree)
void bubble_down_bubbles(BubbleTree& bubble_tree);

// Enumerate Cactus bubbles.  Interface (and output on DAGs)
// identical to superbubbles()
// Note: input graph will be sorted (as done for superbubbles())
map<pair<id_t, id_t>, vector<id_t> > cactusbubbles(VG& graph);

// Convert back from Cactus to VG
// (to, for example, display using vg view)
// todo: also provide mapping info to get nodes embedded in cactus components
VG cactus_to_vg(stCactusGraph* cactus_graph);

// Convert vg into vg formatted cactus representation
// Input graph must be sorted!
VG cactusify(VG& graph);

}


#endif
