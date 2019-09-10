#ifndef VG_CHUNKER_HPP_INCLUDED
#define VG_CHUNKER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "lru_cache.h"
#include "vg.hpp"
#include "xg.hpp"
#include "json2pb.h"
#include "region.hpp"

namespace vg {

using namespace std;


/** Chunk up a graph along a path, using a given number of
 * context expansion steps to fill out the chunks.  Most of the 
 * work done by exising xg functions.
 */
class PathChunker {

public:
   
    // graph used for all path splitting and subgraphing operations
    const PathPositionHandleGraph* graph;

    PathChunker(const PathPositionHandleGraph* graph = NULL);
    ~PathChunker();

    /** Extract subgraph corresponding to given path region into its 
     * own vg graph, and send it to out_stream.  The boundaries of the
     * extracted graph (which can be different because we expand context and don't
     * cut nodes) are written to out_region.  If forward_only set, context
     * is only expanded in the forward direction
     *
     * NOTE: we follow convention of Region coordinates being 0-based 
     * inclusive. 
     * */
    void extract_subgraph(const Region& region, int context, int length, bool forward_only,
                          VG& subgraph, Region& out_region);

    /**
     * Like above, but use (inclusive) id range instead of region on path.
     */
    void extract_id_range(vg::id_t start, vg::id_t end, int context, int length, bool forward_only,
                         VG& subgraph, Region& out_region);

    /**
     * Get a set of all edges in the graph along a path region (to check for discontinuities later on)
     */ 
    set<pair<pair<id_t, bool>, pair<id_t, bool>>> get_path_edge_index(step_handle_t start_step,
                                                                      step_handle_t end_step, int context) const;

};


}

#endif
