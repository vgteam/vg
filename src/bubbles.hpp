#ifndef VG_BUBBLES_HPP_INCLUDED

#define VG_BUBBLES_HPP_INCLUDED

#include <vector>
#include <map>

#include "types.hpp"
#include "utility.hpp"
#include "nodeside.hpp"

extern "C" {
    typedef struct _stCactusGraph stCactusGraph;
    typedef struct _stCactusNode stCactusNode;
    typedef struct _stList stList;
}

using namespace std;

namespace vg {

class VG;

// We use this when talking to Cactus.
struct CactusSide {
    int64_t node;
    bool is_end;
};

// Consolidate bubble finding code here to keep vg class size from getting even
// more out of hand

// ULTRA BUBBLES
// todo : refactor harmonize interface with stuff in deconstruct
// and superbubbles in general
struct Bubble {
    NodeSide start;
    NodeSide end;
    vector<id_t> contents;
    /// cactus now gives us chaining information, stick here for now
    /// so chain_offsets[i]-chain_offsets[i+1] mark the range
    /// of children in chain i.  existing code that doesn't use
    /// chains will be unaffected. 
    vector<int> chain_offsets;
    /// Set to false if the graph contains cycles, and true if it is a DAG.
    bool dag;
    /// Set to true if the graph contains any tips.
    bool tips;
};

typedef Tree<Bubble> BubbleTree;

// Convert VG to Cactus Graph. Takes a list of path names to use to find
// telomeres if present in a connected component.
// Notes:
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns a Cactus graph, and a list of stCactusEdgeEnd* telomeres, in pairs of adjacent items.
pair<stCactusGraph*, stList*> vg_to_cactus(VG& graph, const unordered_set<string>& hint_paths);

// Return the hierchical cactus decomposition
// The root node is meaningless.  Its children are the top level chains.
// The returned tree must be deleted
BubbleTree* ultrabubble_tree(VG& graph);

// Enumerate ultra bubbles.  Interface (and output on DAGs)
// identical to superbubbles()
// Note: input graph will be sorted (as done for superbubbles())
map<pair<id_t, id_t>, vector<id_t> > ultrabubbles(VG& graph);

// Convert back from Cactus to VG
// (to, for example, display using vg view)
// todo: also provide mapping info to get nodes embedded in cactus components
VG cactus_to_vg(stCactusGraph* cactus_graph);

// Convert vg into vg formatted cactus representation
// Input graph must be sorted!
VG cactusify(VG& graph);

}


#endif
