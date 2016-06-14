#ifndef BUBBLES_H

#define BUBBLES_H

#include <vector>
#include <map>

#include "types.hpp"
#include "DetectSuperBubble.hpp"
extern "C" {
#include "sonLib.h"
#include "stCactusGraphs.h"
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


}


#endif
