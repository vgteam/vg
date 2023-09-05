#ifndef VG_CACTUS_HPP_INCLUDED

#define VG_CACTUS_HPP_INCLUDED

/** \file
 * cactus.hpp: Wrapper utility functions for for pinchesAndCacti
 */

#include <vector>
#include <map>

#include "types.hpp"
#include "utility.hpp"
#include "nodeside.hpp"
#include "vg.hpp"

extern "C" {
#include <sonLib/sonLib.h>
#include <sonLib/stCactusGraphs.h>
}

using namespace std;

namespace vg {

// We use this when talking to Cactus.
struct CactusSide {
    int64_t node;
    bool is_end;
};

// Convert VG to Cactus Graph. Takes a list of path names to use to find
// telomeres if present in a connected component.
// If we know the graph is a single weakly connected component, single_component can
// be set to ture to avoid recomputing components.
// Notes:
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns a Cactus graph, and a list of stCactusEdgeEnd* telomeres, in pairs of adjacent items.
pair<stCactusGraph*, stList*> handle_graph_to_cactus(const PathHandleGraph& graph, const unordered_set<string>& hint_paths,
                                                     bool single_component = false);

// Convert back from Cactus to VG
// (to, for example, display using vg view)
// todo: also provide mapping info to get nodes embedded in cactus components
VG cactus_to_vg(stCactusGraph* cactus_graph);

// Convert vg into vg formatted cactus representation
// Input graph must be sorted!
VG cactusify(VG& graph);

}


#endif
