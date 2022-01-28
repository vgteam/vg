/**
 * \file register_loader_saver_gfa.cpp
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gfa.hpp"
#include "algorithms/gfa_to_handle.hpp"
#include "gfa.hpp"

#include "handle.hpp"
#include "bdsg/packed_graph.hpp"
#include "save_handle_graph.hpp"

#include <iostream>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

// derived from Registry::sniff_magic() in libvgio/src/registry.cpp
// note that we have a limited number of characters to work with so
// there's no way to do something extremely precise or general.
// this only works when we can seek around a little bet, which the calling
// logic in libvgio verifies...
static bool sniff_gfa(istream& stream) {
    if (!stream) {
        // Can't read anything, so obviously it can't match.
        return false;
    }
    
    // Work out how many characters to try and sniff.
    // We require that our C++ STL can do this much ungetting, even though the
    // standard guarantees absolutely no ungetting.
    const size_t to_sniff = 8;
    
    // Allocate a buffer
    char buffer[to_sniff];
    // Have a cursor in the buffer
    size_t buffer_used = 0;

    while (stream.peek() != EOF && buffer_used < to_sniff) {
        // Until we fill the buffer or would hit EOF, fill the buffer
        buffer[buffer_used] = (char) stream.get();
        buffer_used++;
    }

    for (size_t i = 0; i < buffer_used; i++) {
        // Now unget all the characters again.
        // C++11 says we can unget from EOF.
        stream.unget();
        if (!stream) {
            // We did something the stream disliked.
            throw runtime_error("Ungetting failed after " + to_string(i) + " characters");
        }
    }
    
    // Now all the characters are back in the stream.
    
    if (!stream) {
        // We reached EOF when sniffing the magic. We managed to unget
        // everything (maybe the file is empty). But we need to clear errors on
        // the stream so it is like it was when we started.
        stream.clear();
    }

    if (buffer_used < 2) {
        // Todo: GFAs can be empty -- may want to return true for empty files?
        // As it stands, we require "H\n" at the very least.
        return false;
    }

    // We are not accepting leading whitespaces.  There is no way to do this generally
    // with our limited buffer.  Also, this check is not bulletproof, so best to sniff gfas
    // as a last resort.

    // Check for a header line
    if (buffer[0] == 'H' &&
        (buffer[1] == '\n' ||
         (buffer[1] == '\t' && buffer_used >= 8 && strncmp(buffer + 2, "VN:Z:", 5) == 0))) {
        return true;
    }

    // Check for any other type of line, looking for a <record type> <TAB> <printable char>
    if ((buffer[0] == 'S' || buffer[0] == 'L' || buffer[0] == 'P' || buffer[0] == 'W') &&
        buffer[1] == '\t' && buffer_used > 3 && isprint(buffer[2])) {
        return true;
    }
    
    return false;
}

void register_loader_saver_gfa() {
    Registry::register_bare_loader_saver_with_header_check<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("GFA", sniff_gfa, [](istream& input, const string& filename) -> void* {
        // Allocate a PackedGraph
        GFAHandleGraph* gfa_graph = new GFAHandleGraph();
        
        try {
        
            if (!filename.empty() && filename != "-") {
                // Load it from a file
                algorithms::gfa_to_path_handle_graph(filename, gfa_graph, &gfa_graph->gfa_id_space);
            } else {
                // Load it from the stream, falling back to temp file if necessary
                algorithms::gfa_to_path_handle_graph_stream(input, gfa_graph, &gfa_graph->gfa_id_space);
            }
            // Make sure the node ID to sequence space translation is ready if anybody wants it.
            gfa_graph->gfa_id_space.invert_translation();
        
        } catch (algorithms::GFAFormatError& e) {
            // There is something wrong with the input GFA file.
            // Explain that without a long stack trace, and bail.
            cerr << "error[register_loader_saver_gfa] GFA ";
            if (!filename.empty() && filename != "-") {
                cerr << "file " << filename;
            } else {
                cerr << "stream";
            }
            cerr << " is corrupt and cannot be loaded." << endl;
            cerr << e.what() << endl;
            exit(1);
        }
        
        // Return it so the caller owns it.
        return (void*) gfa_graph;
    }, [](const void* gfa_graph_void, ostream& output) {
        // Cast to GFAHandleGraph and serialize to the stream.
        assert(gfa_graph_void != nullptr);
        graph_to_gfa((const GFAHandleGraph*)gfa_graph_void, output);
    });
}

}

}

