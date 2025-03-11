/**
 * \file register_loader_saver_gbwtgraph.cpp
 * Defines IO for GBWTGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbwtgraph.hpp"

#include <gbwtgraph/gbwtgraph.h>

#include <arpa/inet.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbwtgraph() {

    // Use the `SerializableHandleGraph` magic number.
    gbwtgraph::GBWTGraph empty;
    std::uint32_t magic_number = htonl(empty.get_magic_number());
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gbwtgraph::GBWTGraph>("GBWTGraph", magic_string, [](istream& input) -> void* {
        gbwtgraph::GBWTGraph* graph = new gbwtgraph::GBWTGraph();

        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Throw an exception or fail silently if reading a Simple-SDS input fails.
        // * Exit with std::exit() or fail silently if reading an SDSL input fails.
        // The exceptions are derived from std::runtime_error.
        graph->deserialize(input);
        
        // Return the graph so the caller owns it.
        return static_cast<void*>(graph);
    }, [](const void* graph_void, ostream& output) {
        assert(graph_void != nullptr);
        // Serialize in the SDSL format, which is larger than the simple-sds format but faster to load.
        // If we want to use the simple-sds format, we can serialize GBZ instead.
        static_cast<const gbwtgraph::GBWTGraph*>(graph_void)->serialize(output);
    });
}

}

}
