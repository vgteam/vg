/**
 * \file register_loader_saver_gbwtgraph.cpp
 * Defines IO for GBWTGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbwtgraph.hpp"

#include <gbwtgraph/gbwtgraph.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbwtgraph() {
    Registry::register_bare_loader_saver<gbwtgraph::GBWTGraph>("GBWTGraph", [](istream& input) -> void* {
        gbwtgraph::GBWTGraph* graph = new gbwtgraph::GBWTGraph();
        graph->deserialize(input);
        
        // Return the graph so the caller owns it.
        return static_cast<void*>(graph);
    }, [](const void* graph_void, ostream& output) {
        assert(graph_void != nullptr);
        static_cast<const gbwtgraph::GBWTGraph*>(graph_void)->serialize(output);
    });
}

}

}
