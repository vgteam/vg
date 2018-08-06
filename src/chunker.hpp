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
#include "gam_index.hpp"

namespace vg {

using namespace std;


/** Chunk up a graph along a path, using a given number of
 * context expansion steps to fill out the chunks.  Most of the 
 * work done by exising xg functions. For gams, the rocksdb
 * index is also required. 
 */
class PathChunker {

public:
   
    // xg index used for all path splitting and subgraphing operations
    xg::XG* xg;
    // number of gams to write at once
    size_t gam_buffer_size = 100;

    PathChunker(xg::XG* xg = NULL);
    ~PathChunker();

    /** Extract subgraph corresponding to given path region into its 
     * own vg graph, and send it to out_stream.  The boundaries of the
     * extracted graph (which can be different because we expand context and don't
     * cut nodes) are written to out_region.  If forward_only set, context
     * is only expanded in the forward direction
     *
     * NOTE: we follow convention of Region coordinates being 1-based 
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
     * Extract all alignments that touch a node in a subgraph and write them 
     * to an output stream using the given GAM index (and this->gam_buffer_size)
     */
    void extract_gam_for_subgraph(VG& subgraph, GAMIndex::cursor_t& cursor, const GAMIndex& index,
                                  ostream* out_stream, bool only_fully_contained = false);                                     
    
    /** 
     * More general interface used by above two functions.
     * Extract GAM data for a series of inclusive contiguous ranges of IDs.
     * The ranges must be sorted and coalesced.
     */
    void extract_gam_for_ids(const vector<pair<vg::id_t, vg::id_t>>& graph_id_ranges, GAMIndex::cursor_t& cursor, const GAMIndex& index,
                             ostream* out_stream, bool only_fully_contained = false);
    
};


}

#endif
