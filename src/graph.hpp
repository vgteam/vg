#ifndef VG_GRAPH_HPP_INCLUDED
#define VG_GRAPH_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "handle.hpp"
#include <set>
#include <algorithm>

namespace vg {

using namespace std;

/// remove duplicates and sort by id
void sort_by_id_dedup_and_clean(Graph& graph);

/// remove duplicate nodes and edges
void remove_duplicates(Graph& graph);

/// remove duplicate edges
void remove_duplicate_edges(Graph& graph);

/// remove duplicate nodes
void remove_duplicate_nodes(Graph& graph);

/// remove edges that link to a node that is not in the graph
void remove_orphan_edges(Graph& graph);

/// order the nodes and edges in the graph by id
void sort_by_id(Graph& graph);

/// order the nodes in the graph by id
void sort_nodes_by_id(Graph& graph);

/// order the edges in the graph by id pairs
void sort_edges_by_id(Graph& graph);

/// returns true if the graph is id-sortable (no reverse links)
bool is_id_sortable(const Graph& graph);

/// returns true if we find an edge that may specify an inversion
bool has_inversion(const Graph& graph);

/// clean up doubly-reversed edges
void flip_doubly_reversed_edges(Graph& graph);

// transfer data from a HandleGraph into an empty Graph
void from_handle_graph(const HandleGraph& from, Graph& to);

// transfer data from a PathHandleGraph into an empty Graph
void from_path_handle_graph(const PathHandleGraph& from, Graph& to);
    
}

#endif
