/**
 * \file register_loader_saver_odgi.cpp
 * Defines IO for a PackedGraph from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_odgi.hpp"

#include "handle.hpp"
#include "bdsg/odgi.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_odgi() {

    // Convert the ODGI SerializableHandleGraph magic number to a string
    bdsg::ODGI empty;
    // Make sure it is in network byte order
    uint32_t new_magic_number = htonl(empty.get_magic_number());
    // Load all 4 characters of it into a string
    string new_magic((char*)&new_magic_number, 4);

    Registry::register_bare_loader_saver_with_magic<bdsg::ODGI, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("ODGI", new_magic, [](istream& input) -> void* {
        // Allocate an ODGI graph
        bdsg::ODGI* odgi = new bdsg::ODGI();
        
        // Load it
        odgi->deserialize(input);
        
        // Return it so the caller owns it.
        return (void*) odgi;
    }, [](const void* odgi_void, ostream& output) {
        // Cast to ODGI and serialize to the stream.
        assert(odgi_void != nullptr);
        ((bdsg::ODGI*) odgi_void)->serialize(output);
    });
}

}

}

