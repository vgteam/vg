#ifndef VG_IO_JSON2GRAPH_HPP_INCLUDED
#define VG_IO_JSON2GRAPH_HPP_INCLUDED

/**
 * \file json2graph.hpp
 * Load a graph from JSON.
 */

#include <handlegraph/algorithms/copy_graph.hpp>

#include <vg/io/json2pb.h>
#include "../vg.hpp"

namespace vg {

namespace io {


/// Load a JSON string into a graph. The string must be a single JSON object.
inline void json2graph(const std::string& json, MutablePathMutableHandleGraph* dest) {
    // Load as a Protobuf message
    Graph g;
    json2pb(g, json);
    
    // Wrap the graph in a HandleGraph
    VG graph(g);

    // And copy to the destination.
    handlegraph::algorithms::copy_path_handle_graph(&graph, dest);
}

}

}

#endif
