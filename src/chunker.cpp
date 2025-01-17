#include <iostream>
#include <unordered_set>
#include <vg/io/stream.hpp>
#include "chunker.hpp"
#include "algorithms/subgraph.hpp"
#include "vg.hpp"
#include "clip.hpp"
#include "crash.hpp"

//#define debug

namespace vg {

using namespace std;

PathChunker::PathChunker(const PathPositionHandleGraph* graph) : graph(graph) {
    
}

PathChunker::~PathChunker() {

}

void PathChunker::extract_subgraph(const Region& region, int64_t context, int64_t length, bool forward_only,
                                   MutablePathMutableHandleGraph& subgraph, Region& out_region) {
    
    // extract our path range into the graph
    // TODO: Handle incoming names with subranges when they aren't exactly the names of paths we have.
    path_handle_t path_handle = graph->get_path_handle(region.seq);
    step_handle_t start_step = graph->get_step_at_position(path_handle, region.start);
    handle_t start_handle = graph->get_handle_of_step(start_step);
    step_handle_t end_step = graph->get_step_at_position(path_handle, region.end);    
    handle_t end_handle = graph->get_handle_of_step(end_step);

    // If the target path was itself a subrange, we need to know the base path range that the requested range of that subrange was.
    subrange_t base_path_subrange = graph->get_subrange(path_handle);

    // Get end-exclusive 0-based subrange we want
    subrange_t target_subrange {region.start, region.end + 1};
    if (base_path_subrange != PathMetadata::NO_SUBRANGE) {
        // Budge it over by the coordinates of what the region is in
        target_subrange.first += base_path_subrange.first;
        target_subrange.second += base_path_subrange.first;
    }
    

#ifdef debug
#pragma omp critical(cerr)
    {
        cerr << "extracting subgraph range for " << region.seq << ":" << region.start << "-" << region.end
             << ", wich maps to handle range " << graph->get_id(start_handle) << ":" << graph->get_is_reverse(start_handle) << "-"
             << graph->get_id(end_handle) << ":" << graph->get_is_reverse(end_handle) << endl;
    }
#endif

    step_handle_t end_plus_one_step = graph->has_next_step(end_step) ? graph->get_next_step(end_step) : graph->path_end(path_handle) ;
    for (step_handle_t step = start_step; step != end_plus_one_step; step = graph->get_next_step(step)) {
        handle_t step_handle = graph->get_handle_of_step(step);
        if (graph->get_is_reverse(step_handle)) {
            step_handle = graph->flip(step_handle);
        }
        if (!subgraph.has_node(graph->get_id(step_handle))) {
            subgraph.create_handle(graph->get_sequence(step_handle), graph->get_id(step_handle));
        }
    };
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    if (context > 0) {
        algorithms::expand_subgraph_by_steps(*graph, subgraph, context, forward_only);
    }
    if (length > 0) {
        algorithms::expand_subgraph_by_length(*graph, subgraph, context, forward_only);
    }
    else if (context == 0 && length == 0) {
        algorithms::add_connecting_edges_to_subgraph(*graph, subgraph);
    }
    algorithms::add_subpaths_to_subgraph(*graph, subgraph); 
   
    // Now we need to figure out how we expanded the target path region we
    // asked for.
    //
    // We don't just want the lowest and highest bounds of any subpath, we want
    // the lowest and highest bound of the subpaths that actually overlap the
    // targeted region.

    // Find the lowest and highest offsets visited by any subpath of the target path we extracted on.
    PathSense sense = graph->get_sense(path_handle);
    std::string sample = graph->get_sample_name(path_handle);
    std::string locus = graph->get_locus_name(path_handle);
    size_t haplotype = graph->get_haplotype(path_handle);
    size_t phase_block = graph->get_phase_block(path_handle);
    

    // Find the outer bounds of selected subpaths of the target path
    size_t min_start = std::numeric_limits<size_t>::max();
    size_t max_end = 0;

    subgraph.for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t subpath) {
        if (subgraph.get_haplotype(subpath) != haplotype || subgraph.get_phase_block(subpath) != phase_block) {
            // Skip this subpath since it's not the right phase/fragment
            return true;
        }

        subrange_t subpath_subrange = subgraph.get_subrange(subpath);
        if (subpath_subrange == PathMetadata::NO_SUBRANGE) {
            // Fill in a 0 start
            subpath_subrange.first = 0;
        }
        if (subpath_subrange.second == PathMetadata::NO_END_POSITION) {
            // Compute a length and use that to get the end.
            // TODO: Sniff for an efficient/available get_path_length.
            size_t path_length = 0;
            for (handle_t handle : subgraph.scan_path(subpath)) {
                path_length += subgraph.get_length(handle);
            }
            subpath_subrange.second = subpath_subrange.first + path_length;
        }

        if (subpath_subrange.first >= target_subrange.second || subpath_subrange.second <= target_subrange.first) {
            // This subpath doesn't actually overlap the selected target base path
            // subrange (which is 0-based, end-exclusive), and so shouldn't count
            // for extending the selected region along the target path.
            return true;
        }

        // Min/max in the subrange bounds
        min_start = std::min(min_start, subpath_subrange.first);
        max_end = std::max(max_end, subpath_subrange.second);

        return true;
    });

    // TODO: We assume we actually found some of the target path
    crash_unless(min_start != std::numeric_limits<size_t>::max());

    // Hackily remove source path subrange offsets if any
    subrange_t source_subrange = graph->get_subrange(path_handle);
    if (source_subrange != PathMetadata::NO_SUBRANGE) {
        // If we requested something on this path region, we can't handle
        // finding part of an earlier path region.
        // TODO: Handle it.
        crash_unless(min_start <= source_subrange.first);
        min_start -= source_subrange.first;
        max_end -= source_subrange.first;
    }

    // We can't represent a region with a 0 end-exclusive coordinate.
    crash_unless(max_end != 0);

    // Produce the output region. Convert coordinates to be 0-based, end-inclusive.
    out_region.seq = region.seq;
    out_region.start = min_start;
    out_region.end = max_end - 1;
}

void PathChunker::extract_snarls(const Region& region, SnarlManager& snarl_manager, MutablePathMutableHandleGraph& subgraph) {

    // copy over the path extraction code from above:

    // extract our path range into the graph
    path_handle_t path_handle = graph->get_path_handle(region.seq);
    step_handle_t start_step = graph->get_step_at_position(path_handle, region.start);
    handle_t start_handle = graph->get_handle_of_step(start_step);
    step_handle_t end_step = graph->get_step_at_position(path_handle, region.end);    
    handle_t end_handle = graph->get_handle_of_step(end_step);

#ifdef debug
#pragma omp critical(cerr)
    {
        cerr << "extracting subgraph range for " << region.seq << ":" << region.start << "-" << region.end
             << ", wich maps to handle range " << graph->get_id(start_handle) << ":" << graph->get_is_reverse(start_handle) << "-"
             << graph->get_id(end_handle) << ":" << graph->get_is_reverse(end_handle) << endl;
    }
#endif

    step_handle_t end_plus_one_step = graph->has_next_step(end_step) ? graph->get_next_step(end_step) : graph->path_end(path_handle) ;
    for (step_handle_t step = start_step; step != end_plus_one_step; step = graph->get_next_step(step)) {
        handle_t step_handle = graph->get_handle_of_step(step);
        if (graph->get_is_reverse(step_handle)) {
            step_handle = graph->flip(step_handle);
        }
        if (!subgraph.has_node(graph->get_id(step_handle))) {
            subgraph.create_handle(graph->get_sequence(step_handle), graph->get_id(step_handle));
        }
    }

    // now fill in the snarls using the vg clip api
    // todo: we can specifiy multiple regions here
    visit_contained_snarls(graph, {region}, snarl_manager, false,
                           [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                               int64_t start_node, int64_t end_node, bool steps_reversed,
                               const Region* containing_region) {

                               pair<unordered_set<id_t>, unordered_set<edge_t> > snarl_contents = snarl_manager.deep_contents(snarl, *graph, true);
                               for (id_t snarl_node : snarl_contents.first) {
                                   if (!subgraph.has_node(snarl_node)) {
                                       subgraph.create_handle(graph->get_sequence(graph->get_handle(snarl_node)), snarl_node);
                                   }
                               }
                               
                           });

    // now fill in the edges
    algorithms::add_connecting_edges_to_subgraph(*graph, subgraph);

    // now fill in the paths
    algorithms::add_subpaths_to_subgraph(*graph, subgraph);
}

void PathChunker::extract_path_component(const string& path_name, MutablePathMutableHandleGraph& subgraph, Region& out_region) {
    unordered_set<nid_t> path_ids;

    path_handle_t path_handle = graph->get_path_handle(path_name);
    for (handle_t handle : graph->scan_path(path_handle)) {
        path_ids.insert(graph->get_id(handle));
    }
    
    extract_component(path_ids, subgraph);
    out_region.seq = path_name;
}

void PathChunker::extract_component(const unordered_set<nid_t>& node_ids, MutablePathMutableHandleGraph& subgraph) {

    for (nid_t node_id : node_ids) {
        subgraph.create_handle(graph->get_sequence(graph->get_handle(node_id)), node_id);
    }

    algorithms::expand_subgraph_by_steps(*graph, subgraph, numeric_limits<uint64_t>::max());
    algorithms::add_subpaths_to_subgraph(*graph, subgraph);
}

void PathChunker::extract_id_range(vg::id_t start, vg::id_t end, int64_t context, int64_t length,
                                   bool forward_only, MutablePathMutableHandleGraph& subgraph,
                                   Region& out_region) {

    for (vg::id_t i = start; i <= end; ++i) {
        subgraph.create_handle(graph->get_sequence(graph->get_handle(i)), i);
    }

    // expand the context and get path information
    // if forward_only true, then we only go forward.
    algorithms::expand_subgraph_by_steps(*graph, subgraph, context, forward_only);
    if (length) {
        algorithms::expand_subgraph_by_length(*graph, subgraph, context, forward_only);
    }
    algorithms::add_subpaths_to_subgraph(*graph, subgraph);

    // build the vg
    out_region.start = subgraph.min_node_id();
    out_region.end = subgraph.max_node_id();
}

set<pair<pair<id_t, bool>, pair<id_t, bool>>> PathChunker::get_path_edge_index(step_handle_t start_step,
                                                                               step_handle_t end_step, int64_t context) const {
    // we don't use handles as we're going to use this structure to compare edges across different graphs
    set<pair<pair<id_t, bool>, pair<id_t, bool>>> path_edges;

    function<void(step_handle_t)> add_edge = [&](step_handle_t step) {
        step_handle_t next = graph->get_next_step(step);
        edge_t edge = graph->edge_handle(graph->get_handle_of_step(step), graph->get_handle_of_step(next));
        path_edges.insert(make_pair(make_pair(graph->get_id(edge.first), graph->get_is_reverse(edge.first)),
                                    make_pair(graph->get_id(edge.second), graph->get_is_reverse(edge.second))));
    };

    // edges from left context
    int i = 0;
    for (step_handle_t step = start_step; graph->has_previous_step(step) && i <= context;
         step = graph->get_previous_step(step), ++i) {
        add_edge(graph->get_previous_step(step));
    }

    // edges from range
    for (step_handle_t step = start_step; step != end_step; step = graph->get_next_step(step)) {
        if (graph->has_next_step(step)) {
            add_edge(step);
        }
    }

    // edges from right context
    i = 0;
    for (step_handle_t step = end_step; graph->has_next_step(step) && i <= context;
         step = graph->get_next_step(step), ++i) {
        add_edge(step);
    }

    return path_edges;
}


}
