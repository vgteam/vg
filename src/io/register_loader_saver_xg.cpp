/**
 * \file register_loader_saver_xg.cpp
 * Defines IO for an XG index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_xg.hpp"

#include "handle.hpp"
#include "xg.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_xg() {
    // Register to load with either old or new SerializableHandleGraph-managed
    // XG magic number sequences, in addition to the tag.
    Registry::register_bare_loader_saver_with_magics<xg::XG, PathPositionHandleGraph, PathHandleGraph, HandleGraph>("XG",
        {"XG", "NGXG"}, [](istream& input) -> void* {
        // Allocate an XG
        xg::XG* index = new xg::XG();
        
        // Load it
        index->deserialize(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to XG and serialize to the stream.
        assert(index_void != nullptr);
        ((const xg::XG*) index_void)->serialize(output);
    });
}

}

}

