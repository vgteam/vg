#include <iostream>
#include <unordered_set>
#include <vg/io/stream.hpp>
#include "chunker.hpp"
#include "algorithms/subgraph.hpp"

//#define debug

namespace vg {

using namespace std;

PathChunker::PathChunker(const PathPositionHandleGraph* graph) : graph(graph) {
    
}

PathChunker::~PathChunker() {

}

void PathChunker::extract_subgraph(const Region& region, int context, int length, bool forward_only,
                                   VG& subgraph, Region& out_region) {
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
        
    // build the vg of the subgraph
    subgraph.remove_orphan_edges();

    // get our range endpoints before context expansion
    list<mapping_t>& mappings = subgraph.paths.get_path(region.seq);
    assert(!mappings.empty());
    size_t mappings_size = mappings.size();
    int64_t input_start_node = graph->get_id(start_handle);
    int64_t input_end_node = graph->get_id(end_handle);

#ifdef debug
#pragma omp critical(cerr)
    {
        cerr << "Path range in expanded subgraph is " << *mappings.begin() << "-" << *mappings.rbegin() << endl;
    }
#endif

    // replaces old xg position_in_path() to check node counts in path
    function<vector<step_handle_t>(const PathHandleGraph&, handle_t, path_handle_t)> path_steps_of_handle =
        [] (const PathHandleGraph& graph, handle_t handle, path_handle_t path_handle) {
        vector<step_handle_t> node_steps = graph.steps_of_handle(handle);
        vector<step_handle_t> node_path_steps;
        for (auto step : node_steps) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                node_path_steps.push_back(step);
            }
        }
        return node_path_steps;
    };
    
    // we have no direct way of getting our steps out of the subgraph, so we
    // go through node ids.  the problem is that cycles can introduce
    // ambiguity.  we check for that here (only to punt on it later)
    vector<step_handle_t> start_node_path_steps = path_steps_of_handle(*graph, start_handle, path_handle);
    vector<step_handle_t> end_node_path_steps = path_steps_of_handle(*graph, end_handle, path_handle);
    bool end_points_on_cycle = start_node_path_steps.size() > 1 || end_node_path_steps.size() > 1;
    
    // keep track of the edges in our original path
    set<pair<pair<id_t, bool>, pair<id_t, bool>>> path_edge_set =
        // walking out with the context length (as supported below) won't always work as expansion
        // can grab an arbitrary amount of path regardless of context.  so we load up the entire path:
        // (todo: could sniff out limits from subgraph...)
        get_path_edge_index(graph->path_begin(path_handle), graph->path_back(path_handle), std::max(context, length));
    
    // the distance between them and the nodes in our input range
    size_t left_padding = 0;
    size_t right_padding = 0;
    // do we need to rewrite back to our graph?
    bool rewrite_paths = false;
    
    if (!end_points_on_cycle) {
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
                handle_t  prev_handle = subgraph.get_handle(prev_it->node_id(),
                                                            prev_it->is_reverse());
                handle_t cur_handle = subgraph.get_handle(cur_it->node_id(),
                                                          cur_it->is_reverse());
                edge_t edge = subgraph.edge_handle(cur_handle, prev_handle);
                if (!path_edge_set.count(make_pair(make_pair(subgraph.get_id(edge.first), subgraph.get_is_reverse(edge.first)),
                                                   make_pair(subgraph.get_id(edge.second), subgraph.get_is_reverse(edge.second))))) {
#ifdef debug
#pragma omp critical(cerr)
                    {
                        cerr << "found discontinuity between when left scanning path in subgraph: " << *cur_it << " and " << *prev_it << endl;

                    }
#endif
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
            handle_t  prev_handle = subgraph.get_handle(prev_it->node_id(),
                                                        prev_it->is_reverse());
            handle_t cur_handle = subgraph.get_handle(cur_it->node_id(),
                                                      cur_it->is_reverse());
            edge_t edge = subgraph.edge_handle(prev_handle, cur_handle);
            if (!path_edge_set.count(make_pair(make_pair(subgraph.get_id(edge.first), subgraph.get_is_reverse(edge.first)),
                                               make_pair(subgraph.get_id(edge.second), subgraph.get_is_reverse(edge.second))))) {
#ifdef debug
#pragma omp critical(cerr)
                    {
                        cerr << "found discontinuity between when right scanning path in subgraph: " << *prev_it << " and " << *cur_it << endl;

                    }
#endif
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
        for (step_handle_t step = start_step; step != end_plus_one_step; step = graph->get_next_step(step)) {
            handle_t step_handle = graph->get_handle_of_step(step);
            mapping_t mapping;
            mapping.set_node_id(graph->get_id(step_handle));
            mapping.set_is_reverse(graph->get_is_reverse(step_handle));
            mappings.push_back(mapping);
        }
        rewrite_paths = true;
    }

    // Cut our graph so that our reference path end points are graph tips.  This will let the
    // snarl finder use the path to find telomeres.
    path_handle_t sg_path_handle = subgraph.get_path_handle(region.seq);
    Node* start_node = subgraph.get_node(mappings.begin()->node_id());
    auto sg_start_steps = path_steps_of_handle(subgraph, subgraph.get_handle(start_node->id()), sg_path_handle); 
    if (rewrite_paths && sg_start_steps.size() == 1) {
        if (!mappings.begin()->is_reverse() && subgraph.start_degree(start_node) != 0) {
            for (auto edge : subgraph.edges_to(start_node)) {
#ifdef debug
#pragma omp crticial(cerr)
                {
                    cerr << "clipping out edge " << pb2json(*edge) << " in order to make path start a tip" << endl;
                }
#endif
                subgraph.destroy_edge(edge);
            }
        } else if (mappings.begin()->is_reverse() && subgraph.end_degree(start_node) != 0) {
            for (auto edge : subgraph.edges_from(start_node)) {
#ifdef debug
#pragma omp crticial(cerr)
                {
                    cerr << "clipping out edge " << pb2json(*edge) << " in order to make path start a tip" << endl;
                }
#endif
                subgraph.destroy_edge(edge);
            }
        }
    }
    Node* end_node = subgraph.get_node(mappings.rbegin()->node_id());
    auto sg_end_steps = path_steps_of_handle(subgraph, subgraph.get_handle(end_node->id()), sg_path_handle); 
    if (rewrite_paths && sg_end_steps.size() == 1) {
        if (!mappings.rbegin()->is_reverse() && subgraph.end_degree(end_node) != 0) {
            for (auto edge : subgraph.edges_from(end_node)) {
#ifdef debug
#pragma omp crticial(cerr)
                {
                    cerr << "clipping out edge " << pb2json(*edge) << " in order to make path end a tip" << endl;
                }
#endif
                subgraph.destroy_edge(edge);
            }
        } else if (mappings.rbegin()->is_reverse() && subgraph.start_degree(end_node) != 0) {
            for (auto edge : subgraph.edges_to(end_node)) {
#ifdef debug
#pragma omp crticial(cerr)
                {
                    cerr << "clipping out edge " << pb2json(*edge) << " in order to make path end a tip" << endl;
                }
#endif
                subgraph.destroy_edge(edge);
            }
        }
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
    int64_t input_start_pos = graph->get_position_of_step(start_step);
    int64_t input_end_pos = graph->get_position_of_step(end_step);
    out_region.seq = region.seq;
    out_region.start = input_start_pos - left_padding;
    out_region.end = input_end_pos + graph->get_length(end_handle) + right_padding - 1;
}

void PathChunker::extract_id_range(vg::id_t start, vg::id_t end, int context, int length,
                                   bool forward_only, VG& subgraph,
                                   Region& out_region) {

    Graph g;

    for (vg::id_t i = start; i <= end; ++i) {
        Node node;
        node.set_id(i);
        node.set_sequence(graph->get_sequence(graph->get_handle(i)));
        *g.add_node() = node;
    }

    VG vg_g(g);
    
    // expand the context and get path information
    // if forward_only true, then we only go forward.
    algorithms::expand_subgraph_by_steps(*graph, vg_g, context, forward_only);
    if (length) {
        algorithms::expand_subgraph_by_length(*graph, vg_g, context, forward_only);
    }
    algorithms::add_subpaths_to_subgraph(*graph, vg_g);

    // build the vg
    subgraph.extend(vg_g);
    subgraph.remove_orphan_edges();

    out_region.start = subgraph.min_node_id();
    out_region.end = subgraph.max_node_id();
}

set<pair<pair<id_t, bool>, pair<id_t, bool>>> PathChunker::get_path_edge_index(step_handle_t start_step,
                                                                               step_handle_t end_step, int context) const {
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
