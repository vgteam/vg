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
#include "index.hpp"

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
    size_t gam_buffer_size = 1000;

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
    void extract_subgraph(const Region& region, int context, bool forward_only,
                          VG& subgraph, Region& out_region);

    /**
     * Like above, but use (inclusive) id range instead of region on path.
     */
    void extract_id_range(vg::id_t start, vg::id_t end, int context, bool forward_only,
                         VG& subgraph, Region& out_region);

    /** Extract all alignments that touch a node in a subgraph and write them 
     * to an output stream using the rocksdb index (and this->gam_buffer_size) */
    int64_t extract_gam_for_subgraph(VG& subgraph, Index& index, ostream* out_stream,
                                     bool only_fully_contained = false);
    
    /** Like above, but for (inclusive) node id range */
    int64_t extract_gam_for_id_range(vg::id_t start, vg::id_t end, Index& index, ostream* out_stream,
                                     bool only_fully_contained = false);

    /** More general interface used by above two functions */
    int64_t extract_gam_for_ids(const vector<vg::id_t>& graph_ids, Index& index, ostream* out_stream,
                                bool contiguous = false, bool only_fully_contained = false);
    
};


}

#endif
