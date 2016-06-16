#ifndef BUBBLES_H

#define BUBBLES_H

#include <vector>
#include <map>

#include "types.hpp"
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

// Recursive cactus structure: bubbles nested within other bubbles
struct BubbleNode {
    Bubble bubble;
    vector<BubbleNode*> children;
    ~BubbleNode() { for (auto c : children) { delete c; } }
    void for_each_preorder(function<void(Bubble&)> lambda) {
        lambda(bubble);
        for (auto c : children) { c->for_each_preorder(lambda); } }
};
struct BubbleTree {
    BubbleNode* root;
    BubbleTree(BubbleNode* r = 0) : root(r) {}
    ~BubbleTree() { delete root; }
    void for_each_preorder(function<void(Bubble&)> lambda) {
        if (root) root->for_each_preorder(lambda); }
};

// Convert VG to Cactus Graph
// Notes:
//  - source and sink of vg graph obtained by min/max node id
//    so graph should be sorted
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns "root" node as well as graph
pair<stCactusGraph*, stCactusNode*> vg_to_cactus(VG& graph);

// Return the hierchical cactus decomposition
BubbleTree cactusbubble_tree(VG& graph);

// Enumerate Cactus bubbles.  Interface (and output on DAGs)
// identical to superbubbles()
map<pair<id_t, id_t>, vector<id_t> > cactusbubbles(VG& graph);


}


#endif
