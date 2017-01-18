#include <iostream>
#include "stream.hpp"
#include "chunker.hpp"


namespace vg {

using namespace std;
using namespace xg;



PathChunker::PathChunker(xg::XG* xindex) : xg(xindex) {
    
}

PathChunker::~PathChunker() {

}

int64_t PathChunker::extract_subgraph(const Region& region, int context, VG& subgraph) {

    Graph g;

    // convert to 0-based inclusive
    int64_t start = region.start - 1;

    // extract our path range into the graph

    // Commenting out till I can be sure it's not doing weird things to paths
    //xg->get_path_range(region.seq, region.start, region.end - 1, g);

    xg->for_path_range(region.seq, start, region.end - 1, [&](int64_t id) {
            *g.add_node() = xg->node(id);
        });
    
    // expand the context and get path information
    xg->expand_context(g, context, true);
        
    // build the vg
    subgraph.extend(g);
    subgraph.remove_orphan_edges();

    // what node contains our input starting position?
    int64_t input_start_node = xg->node_at_path_position(region.seq, start);

    // start could fall inside a node.  we find out where in the path the
    // 0-offset point of the node is. 
    int64_t input_start_pos = xg->node_start_at_path_position(region.seq, start);
    assert(input_start_pos <= start &&
           input_start_pos + xg->node_length(input_start_node) > start);
    
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
    
    return chunk_start_pos;
}

int64_t PathChunker::extract_gam_for_subgraph(VG& subgraph, Index& index, ostream* out_stream) {

    // Build the set of all the node IDs to operate on
    vector<vg::id_t> graph_ids;
    subgraph.for_each_node([&](Node* node) {
        // Put all the ids in the set
        graph_ids.push_back(node->id());
    });

    // Load all the reads matching the graph into memory
    vector<Alignment> gam_buffer;
    int64_t gam_count = 0;

    function<Alignment&(uint64_t)> write_buffer_elem = [&gam_buffer](uint64_t i) -> Alignment& {
        return gam_buffer[i];
    };

    function<void(const Alignment&)> write_alignment = [&](const Alignment& alignment) {
        // flush our buffer if it's too big
        if (gam_buffer.size() > gam_buffer_size) {
            stream::write(*out_stream, gam_buffer.size(), write_buffer_elem);
            gam_count += gam_buffer.size();
            gam_buffer.clear();
        }
        
        // add to buffer
        gam_buffer.push_back(alignment);
    };

    index.for_alignment_to_nodes(graph_ids, write_alignment);

    // flush buffer 
    if (gam_buffer.size() > 0) {
        stream::write(*out_stream, gam_buffer.size(), write_buffer_elem);
        gam_count += gam_buffer.size();
        gam_buffer.clear();
    }
    
    return gam_count;
}
                
    

}
