#ifndef VG_ALGORITHMS_COPY_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_COPY_GRAPH_HPP_INCLUDED

/**
 * \file copy_graph.hpp
 *
 * Defines algorithms for copying data between handle graphs
 */

#include <iostream>

#include "../handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

    /// Copies the nodes and edges from one graph into another.
    void copy_handle_graph(const HandleGraph* from, MutableHandleGraph* into);
    
    /// Copies the nodes, edges, and paths from one graph into another.
    void copy_path_handle_graph(const PathHandleGraph* from, MutablePathMutableHandleGraph* into);
    
    /// Copies a path from one graph to another. Nodes and edges to support
    /// the path must already exist.
    void copy_path(const PathHandleGraph* from, const path_handle_t& path,
                   MutablePathHandleGraph* into);

}
}

#endif
