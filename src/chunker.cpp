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

void PathChunker::extract_subgraph(const Region& region, int context, bool forward_only,
                                   VG& subgraph, Region& out_region) {

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
    // if forward_only true, then we only go forward.
    xg->expand_context(g, context, true, true, true, !forward_only);
        
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

    out_region.seq = region.seq;
    out_region.start = 1 + chunk_start_pos;
    out_region.end = out_region.start - 1;
    // Is there a better way to get path length? 
    Path output_path = subgraph.paths.path(out_region.seq);
    for (size_t j = 0; j < output_path.mapping_size(); ++j) {
      int64_t op_node = output_path.mapping(j).position().node_id();
      out_region.end += subgraph.get_node(op_node)->sequence().length();
    }
}

void PathChunker::extract_id_range(vg::id_t start, vg::id_t end, int context,
                                   bool forward_only, VG& subgraph, Region& out_region) {

    Graph g;

    for (vg::id_t i = start; i <= end; ++i) {
        *g.add_node() = xg->node(i);
    }
    
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    xg->expand_context(g, context, true, true, true, !forward_only);
        
    // build the vg
    subgraph.extend(g);
    subgraph.remove_orphan_edges();

    out_region.start = subgraph.min_node_id();
    out_region.end = subgraph.max_node_id();
}

int64_t PathChunker::extract_gam_for_subgraph(VG& subgraph, Index& index,
                                              ostream* out_stream,
                                              bool only_fully_contained) {

    // Build the set of all the node IDs to operate on
    vector<vg::id_t> graph_ids;
    subgraph.for_each_node([&](Node* node) {
        // Put all the ids in the set
        graph_ids.push_back(node->id());
    });

    return extract_gam_for_ids(graph_ids, index, out_stream, false,
                               only_fully_contained);
}

int64_t PathChunker::extract_gam_for_id_range(vg::id_t start, vg::id_t end, Index& index,
                                              ostream* out_stream,
                                              bool only_fully_contained) {
    
    vector<vg::id_t> graph_ids;
    for (vg::id_t i = start; i <= end; ++i) {
        graph_ids.push_back(i);
    }
    
    return extract_gam_for_ids(graph_ids, index, out_stream, true,
                               only_fully_contained);
}

int64_t PathChunker::extract_gam_for_ids(const vector<vg::id_t>& graph_ids,
                                         Index& index, ostream* out_stream,
                                         bool contiguous,
                                         bool only_fully_contained) {
  
    // Load all the reads matching the graph into memory
    vector<Alignment> gam_buffer;
    int64_t gam_count = 0;

    function<Alignment&(uint64_t)> write_buffer_elem = [&gam_buffer](uint64_t i) -> Alignment& {
        return gam_buffer[i];
    };

    // We filter out alignments with no nodes in our id range as post-processing
    // until for_alignment_to_nodes fixed to not return such things.
    function<bool(vg::id_t)> check_id;
    unordered_set<vg::id_t>* id_lookup = NULL;
    if (contiguous) {
        check_id = [&](vg::id_t node_id) {
            return node_id >= graph_ids[0] && node_id <= graph_ids[graph_ids.size() - 1];
        };
    } else {
        id_lookup = new unordered_set<vg::id_t>(graph_ids.begin(), graph_ids.end());
        check_id = [&](vg::id_t node_id) {
            return id_lookup->count(node_id) == 1;
        };
    }
    int filter_count = 0;
    function<bool(const Alignment&)> in_range = [&](const Alignment& alignment) {
        if (alignment.has_path()) {
            for (size_t i = 0; i < alignment.path().mapping_size(); ++i) {
                bool check = check_id(alignment.path().mapping(i).position().node_id());
                if (only_fully_contained && check == false) {
                    return false;
                } else if (!only_fully_contained && check == true) {
                    return true;
                }
            }
        }
        if (only_fully_contained) {
            return true;
        } else {
            ++filter_count;
            return false;
        }
    };

    function<void(const Alignment&)> write_alignment = [&](const Alignment& alignment) {
        // flush our buffer if it's too big
        if (gam_buffer.size() > gam_buffer_size) {
            stream::write(*out_stream, gam_buffer.size(), write_buffer_elem);
            gam_count += gam_buffer.size();
            gam_buffer.clear();
        }
        
        // add to buffer
        if (in_range(alignment)) {
            gam_buffer.push_back(alignment);
        }
    };

    index.for_alignment_to_nodes(graph_ids, write_alignment);

    // flush buffer 
    if (gam_buffer.size() > 0) {
        stream::write(*out_stream, gam_buffer.size(), write_buffer_elem);
        gam_count += gam_buffer.size();
        gam_buffer.clear();
    }

    delete id_lookup;
    if (filter_count > 0) {
        cerr << "[vg chunk] Filtered " << filter_count << " erroneous hits from Rocksdb query" << endl;
    }
    
    return gam_count;

}

}
