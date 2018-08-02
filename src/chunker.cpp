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

    out_region.seq = region.seq;
    out_region.start = chunk_start_pos;
    out_region.end = out_region.start - 1;
    // Is there a better way to get path length? 
    Path output_path = subgraph.paths.path(out_region.seq);
    for (size_t j = 0; j < output_path.mapping_size(); ++j) {
      int64_t op_node = output_path.mapping(j).position().node_id();
      out_region.end += subgraph.get_node(op_node)->sequence().length();
    }
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

void PathChunker::extract_gam_for_subgraph(VG& subgraph, GAMIndex::cursor_t& cursor, const GAMIndex& index, 
                                           ostream* out_stream, bool only_fully_contained) {

    // Build the list of of all the node IDs to operate on
    bool contiguous = true;
    vector<vg::id_t> graph_ids;
    subgraph.for_each_node([&](Node* node) {
        // Put all the ids in the list
        graph_ids.push_back(node->id());
    });
    
    // Sort the graph IDs
    std::sort(graph_ids.begin(), graph_ids.end());
    
    // Coalesce them into ranges
    vector<pair<vg::id_t, vg::id_t>> ranges;
    for (auto& id : graph_ids) {
        if (ranges.empty() || ranges.back().second + 1 != id) {
            // We can't glom on to the previous range, so start a new one of just us
            ranges.emplace_back(id, id);
        } else {
            // Extend the previous range to us
            ranges.back().second = id;
        }
    }

    extract_gam_for_ids(ranges, cursor, index, out_stream, only_fully_contained);
}

void PathChunker::extract_gam_for_ids(const vector<pair<vg::id_t, vg::id_t>>& graph_id_ranges,
                                      GAMIndex::cursor_t& cursor, const GAMIndex& index,
                                      ostream* out_stream,  bool only_fully_contained) {


    // Make an output emitter
    stream::ProtobufEmitter<Alignment> emitter(*out_stream);
    
    // Send all the ranges to the index in a combined query.
    // The index exploits their sorted-ness so that it can pull stuff out in at most one pass through each relevant group.
    index.find(cursor, graph_id_ranges, [&](const Alignment& aln) {
        // Emit all the matching alignments.
        // TODO: Set up so we can use the move operation the cursors support
        // Not easy because of https://stackoverflow.com/a/30394755
        emitter.write_copy(aln);
    }, only_fully_contained); // Make sure to restrict to only contained alignments if requested.
}

}
