/**
 * \file register_loader_saver_packed_graph.cpp
 * Defines IO for a PackedGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_packed_graph.hpp"

#include "handle.hpp"
#include "bdsg/packed_graph.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_packed_graph() {
  Registry::register_bare_loader_saver<bdsg::PackedGraph, PathHandleGraph, HandleGraph>("PackedGraph", [](istream& input) -> void* {
        // Allocate a PackedGraph
         bdsg::PackedGraph* packed_graph = new bdsg::PackedGraph();
        
        // Load it
        packed_graph->deserialize(input);
        
        // Return it so the caller owns it.
        return (void*) packed_graph;
    }, [](const void* packed_graph_void, ostream& output) {
        // Cast to PackedGraph and serialize to the stream.
        assert(packed_graph_void != nullptr);
        ((const bdsg::PackedGraph*) packed_graph_void)->serialize(output);
    });
}

}

}

