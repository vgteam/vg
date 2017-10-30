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

// Convert VG to Cactus Graph. Takes a list of path names to use to find
// telomeres if present in a connected component.
// Notes:
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns a Cactus graph, and a list of stCactusEdgeEnd* telomeres, in pairs of adjacent items.
pair<stCactusGraph*, stList*> vg_to_cactus(VG& graph, const unordered_set<string>& hint_paths);

// Convert back from Cactus to VG
// (to, for example, display using vg view)
// todo: also provide mapping info to get nodes embedded in cactus components
VG cactus_to_vg(stCactusGraph* cactus_graph);

// Convert vg into vg formatted cactus representation
// Input graph must be sorted!
VG cactusify(VG& graph);

}


#endif
