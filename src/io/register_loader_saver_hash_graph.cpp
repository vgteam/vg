/**
 * \file register_loader_saver_hash_graph.cpp
 * Defines IO for a HashGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_hash_graph.hpp"

#include "handle.hpp"
#include "bdsg/hash_graph.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_hash_graph() {
  Registry::register_bare_loader_saver<bdsg::HashGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("HashGraph", [](istream& input) -> void* {
        // Allocate a HashGraph
         bdsg::HashGraph* hash_graph = new bdsg::HashGraph();
        
        // Load it
        hash_graph->deserialize(input);
        
        // Return it so the caller owns it.
        return (void*) hash_graph;
    }, [](const void* hash_graph_void, ostream& output) {
        // Cast to HashGraph and serialize to the stream.
        assert(hash_graph_void != nullptr);
        ((const bdsg::HashGraph*) hash_graph_void)->serialize(output);
    });
}

}

}

