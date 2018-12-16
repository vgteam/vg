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

    // extract our path range into the graph
    Graph g;
    xg->for_path_range(region.seq, region.start, region.end, [&](int64_t id, bool) {
            *g.add_node() = xg->node(id);
        });
    
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    xg->expand_context(g, context, true, true, true, !forward_only);
    if (length) {
        xg->expand_context(g, context, true, false, true, !forward_only);
    }
        
    // build the vg of the subgraph
    subgraph.extend(g);
    subgraph.remove_orphan_edges();

    // get our range endpoints before context expansion
    list<mapping_t>& mappings = subgraph.paths.get_path(region.seq);
    size_t mappings_size = mappings.size();
    int64_t input_start_node = xg->node_at_path_position(region.seq, region.start);
    vector<size_t> first_positions = xg->position_in_path(input_start_node, region.seq);
    int64_t input_end_node = xg->node_at_path_position(region.seq, region.end);
    vector<size_t> last_positions = xg->position_in_path(input_end_node, region.seq);

    // the distance between then and the nodes in our input range
    size_t left_padding = 0;
    size_t right_padding = 0;
    // do we need to rewrite back to our graph?
    bool rewrite_paths = false;
    
    // Endpoints not in cycles: we can get our output region directly from xg lookups
    if (first_positions.size() == 1  && last_positions.size() == 1) {
        // start and end of our expanded chunk
        auto start_it = mappings.begin();
        auto end_it = --mappings.end();

        // find our input range in the expanded path. we know these nodes only appear once.
        for (; start_it != mappings.end() && start_it->node_id() != input_start_node; ++start_it);
        for (; end_it != mappings.begin() && end_it->node_id() != input_end_node; --end_it);

        // walk back our start point as we can without rank discontinuities. doesn't matter
        // if we encounter cycles here, because we keep a running path length
        auto cur_it = start_it;
        auto prev_it = cur_it;
        if (prev_it != mappings.begin()) {
            for (; prev_it != mappings.begin(); --prev_it) {
                cur_it = prev_it;
                --cur_it;
                if ((prev_it->rank > 0 || cur_it->rank > 0) && cur_it->rank + 1 != prev_it->rank) {
                    break;
                }
                left_padding += cur_it->length;
            }
        }
        start_it = prev_it;
        // walk forward the end point
        cur_it = end_it;
        prev_it = cur_it;
        for (++cur_it; cur_it != mappings.end(); ++prev_it, ++cur_it) {
            if ((prev_it->rank > 0 || cur_it->rank > 0) && prev_it->rank + 1 != cur_it->rank) {
                break;
            }
            right_padding += cur_it->length;
        }
        end_it = prev_it;

        rewrite_paths = start_it != mappings.begin() || end_it != --mappings.end();
        
        // cut out nodes before and after discontinuity
        mappings.erase(mappings.begin(), start_it);
        mappings.erase(++end_it, mappings.end());
    }
    // We're clipping at a cycle in the reference path.  Just preserve the path as-is from the
    // input region.  
    else {
        mappings.clear();
        xg->for_path_range(region.seq, region.start, region.end, [&](int64_t id, bool rev) {
                mapping_t mapping;
                mapping.set_node_id(id);
                mapping.set_is_reverse(rev);
                mappings.push_back(mapping);
            });
        rewrite_paths = true;
    }
    
    // Sync our updated paths lists back into the Graph protobuf
    if (rewrite_paths) {
        subgraph.paths.rebuild_node_mapping();
        subgraph.paths.rebuild_mapping_aux();
        subgraph.graph.clear_path();
        subgraph.paths.to_graph(subgraph.graph);
    }

    // start could fall inside a node.  we find out where in the path the
    // 0-offset point of the node is. 
    int64_t input_start_pos = xg->node_start_at_path_position(region.seq, region.start);
    int64_t input_end_pos = xg->node_start_at_path_position(region.seq, region.end);
    out_region.seq = region.seq;
    out_region.start = input_start_pos - left_padding;
    out_region.end = input_end_pos + xg->node_length(input_end_node) + right_padding - 1;
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
