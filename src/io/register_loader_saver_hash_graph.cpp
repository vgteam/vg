/**
 * \file register_loader_saver_hash_graph.cpp
 * Defines IO for a HashGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_hash_graph.hpp"
#include "load_proto_to_graph.hpp"

#include "../handle.hpp"
#include <bdsg/hash_graph.hpp>


namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_hash_graph() {
    Registry::register_bare_loader_saver<bdsg::HashGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>(
        "HashGraph", [](istream& input) -> void* {
        
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
    
    // Also register to be able to load Protobuf, by converting to a hash graph on input, if vg::VG is not required.
    // The default implementation for a VG loaded from a file as a handle graph will now be a HashGraph.
    Registry::register_loader<bdsg::HashGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>(
        vector<string>{"VG", ""},
        [](const message_sender_function_t& for_each_message) -> void* {
    
        // Allocate a HashGraph
        bdsg::HashGraph* hash_graph = new bdsg::HashGraph();
        
        // Load into it
        load_proto_to_graph(for_each_message, hash_graph);
        
        // Return it so the caller owns it.
        return (void*) hash_graph;
    });
}

}

}

