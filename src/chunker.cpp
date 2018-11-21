#include <iostream>
#include <unordered_set>
#include "stream.hpp"
#include "chunker.hpp"


namespace vg {

using namespace std;
using namespace xg;



PathChunker::PathChunker(xg::XG* xindex) : xg(xindex) {
    
}

PathChunker::~PathChunker() {

}

void PathChunker::extract_subgraph(const Region& region, int context, int length,
                                   bool forward_only, VG& subgraph, Region& out_region) {

    Graph g;

    // convert to 0-based inclusive
    int64_t start = region.start;

    // extract our path range into the graph

    // Commenting out till I can be sure it's not doing weird things to paths
    //xg->get_path_range(region.seq, region.start, region.end - 1, g);

    xg->for_path_range(region.seq, region.start, region.end, [&](int64_t id) {
            *g.add_node() = xg->node(id);
        });
    
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    xg->expand_context(g, context, true, true, true, !forward_only);
    if (length) {
        xg->expand_context(g, context, true, false, true, !forward_only);
    }
        
    // build the vg
    subgraph.extend(g);
    subgraph.remove_orphan_edges();

    // what node contains our input starting position?
    int64_t input_start_node = xg->node_at_path_position(region.seq, region.start);

    // start could fall inside a node.  we find out where in the path the
    // 0-offset point of the node is. 
    int64_t input_start_pos = xg->node_start_at_path_position(region.seq, region.start);
    assert(input_start_pos <= region.start &&
           input_start_pos + xg->node_length(input_start_node) > region.start);
    
    // find out the start position of the first node in the path in the
    // subgraph.  take the last occurance before the input_start_pos
    // todo: there are probably some cases involving cycles where this breaks
    Path path = subgraph.paths.path(region.seq);
    int64_t chunk_start_node = path.mapping(0).position().node_id();    
    int64_t chunk_start_pos = -1;
    int64_t best_delta = numeric_limits<int64_t>::max();
    vector<size_t> first_positions = xg->position_in_path(chunk_start_node, region.seq);
    for (auto fp : first_positions) {
        int64_t delta = input_start_pos - (int64_t)fp;
        if (delta >= 0 && delta < best_delta) {
            best_delta = delta;
            chunk_start_pos = fp;
        }
    }
    assert(chunk_start_pos >= 0);

    // detect discontinuous reference path in expanded graph (can happen from partially covered deletion)
    list<mapping_t>& mappings = subgraph.paths.get_path(region.seq);
    size_t mappings_size = mappings.size();
    size_t cur_pos = chunk_start_pos;
    size_t out_path_len = mappings.begin()->length;
    auto first_it = mappings.begin();
    auto last_it = mappings.end();
    auto cur_it = mappings.begin();
    ++cur_it;
    for (auto prev_it = mappings.begin(); cur_it != mappings.end(); ++prev_it, ++cur_it) {
        if (prev_it->rank > 0 && cur_it->rank > 0 && prev_it->rank + 1 != cur_it->rank) {
            if (cur_pos < region.start) {
                first_it = cur_it;
                chunk_start_pos = cur_pos;
                out_path_len = 0;
            }  else if (cur_pos > region.end) {
                last_it = prev_it;
                break;
            } else {
                cerr << "Warning: attempting to chunk path " << region.seq << "with discontinuous ranks: "
                     << pb2json(prev_it->to_mapping()) << " -> " << pb2json(cur_it->to_mapping()) << endl;
            }
        }
        cur_pos += prev_it->length;
        out_path_len += cur_it->length;
    }
    // cut out nodes before and after discontinuity
    mappings.erase(mappings.begin(), first_it);
    mappings.erase(last_it, mappings.end());

    // save the paths back to the graph if we changed anything
    if (mappings.size() != mappings_size) {
        subgraph.paths.rebuild_node_mapping();
        subgraph.paths.rebuild_mapping_aux();
        subgraph.graph.clear_path();
        subgraph.paths.to_graph(subgraph.graph);
    }
        
    out_region.seq = region.seq;
    out_region.start = chunk_start_pos;
    out_region.end = chunk_start_pos + out_path_len - 1;
}

void PathChunker::extract_id_range(vg::id_t start, vg::id_t end, int context, int length,
                                   bool forward_only, VG& subgraph, Region& out_region) {

    Graph g;

    for (vg::id_t i = start; i <= end; ++i) {
        *g.add_node() = xg->node(i);
    }
    
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    xg->expand_context(g, context, true, true, true, !forward_only);
    if (length) {
        xg->expand_context(g, context, true, false, true, !forward_only);
    }

    // build the vg
    subgraph.extend(g);
    subgraph.remove_orphan_edges();

    out_region.start = subgraph.min_node_id();
    out_region.end = subgraph.max_node_id();
}

}
