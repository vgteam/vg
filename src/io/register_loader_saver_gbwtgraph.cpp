/**
 * \file register_loader_saver_gbwtgraph.cpp
 * Defines IO for a minimizer index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbwtgraph.hpp"

#include "../gbwt_helper.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbwtgraph() {
    Registry::register_bare_loader_saver<GBWTGraph>("GBWTGraph", [](istream& input) -> void* {
        GBWTGraph* graph = new GBWTGraph();
        graph->deserialize(input);
        
        // Return the graph so the caller owns it.
        return static_cast<void*>(graph);
    }, [](const void* graph_void, ostream& output) {
        assert(graph_void != nullptr);
        static_cast<const GBWTGraph*>(graph_void)->serialize(output);
    });
}

}

}
