/**
 * \file register_loader_saver_vg.cpp
 * Defines IO for a VG graph from stream files of Graph objects.
 */

#include <vg/io/registry.hpp>
#include <vg/io/protobuf_iterator.hpp>
#include "register_loader_saver_vg.hpp"

#include "../vg.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_vg() {
    // We register for "" so we can handle untagged old-style vg files and make them into HandleGraphs
    Registry::register_loader_saver<VG>(vector<string>{"VG", ""},
        [](const message_sender_function_t& for_each_message) -> void* {
        // We have a bit of a control problem.
        // The source function wants to drive; we give it a function of strings, and it calls it with all the strings in turn.
        
        // But the VG also wants to drive; we give it a function to fill Graph objects, and it calls it until it runs out.
        // So we use a new constructor of the VG that we get to drive.
    
        // Allocate a VG and have it call a callback to request all graph chunks.
        VG* vg_graph = new VG([&](const function<void(Graph&)>& consume_graph) {
            // Call the source function with a function that deserializes each message and feeds it to the graph.
            for_each_message([&](const string& serialized_graph) {
                // Parse the message to a Graph
                Graph g;
                if (!ProtobufIterator<Graph>::parse_from_string(g, serialized_graph)) {
                    // Handle bad graphs.
                    // TODO: make this an exception if we ever want to be allowed to continue from this.
                    cerr << "error[register_loader_saver_vg]: invalid Graph message" << endl;
                    exit(1);
                }
                
                // Feed the Graph to the VG.
                // TODO: Is all this callback-and-callforth better than just a loop somewhere?
                consume_graph(g);
            });
        });
        
        // Return it so the caller owns it.
        return (void*) vg_graph;
    }, [](const void* vg_void, const message_consumer_function_t& send_message) {
         assert(vg_void != nullptr);
        
        // Cast to a VG
        // TODO: VG has to do some non-const syncing before it can serialize, so we hackily make a non-const pointer!
        // Fix this when VG learns to serialize itself without modification!
        VG& vg_graph = *((VG*) vg_void);
        
        // Chunk the VG to graphs and spit them out as strings using the given send_message function.
        vg_graph.serialize_to_function([&](const Graph& chunk) {
            // Make the Graph into a string
            string s;
            chunk.SerializeToString(&s);
            
            // Ship it.
            // TODO: make this more move-y.
            send_message(s);
        });
    });
}

}

}

