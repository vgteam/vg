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
};

typedef Tree<Bubble> BubbleTree;

// Convert VG to Cactus Graph
// Notes:
//  - source and sink of vg graph obtained by min/max node id
//    so graph should be sorted
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns "root" node as well as graph
pair<stCactusGraph*, stCactusNode*> vg_to_cactus(VG& graph, id_t source_id, id_t sink_id);

// Get source and sink nodes
pair<id_t, id_t> get_cactus_source_sink(VG& graph);

// Return the hierchical cactus decomposition
BubbleTree cactusbubble_tree(VG& graph);

// By default, bubble X's contents array doesn't have all the nodes
// in its children's contents (ie each node only stored in lowest bubble
// that contains it).  This function will store each node in every
// bubble that contains it (like superbubbles)
// Note: contents wont be in any particular order after
void bubble_up_bubbles(BubbleTree& bubble_tree);

// Enumerate Cactus bubbles.  Interface (and output on DAGs)
// identical to superbubbles()
map<pair<id_t, id_t>, vector<id_t> > cactusbubbles(VG& graph);


}


#endif
