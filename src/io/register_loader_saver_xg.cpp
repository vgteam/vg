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
  Registry::register_bare_loader_saver<xg::XG, PathPositionHandleGraph, PathHandleGraph, HandleGraph>("XG", [](istream& input) -> void* {
        // Allocate an XG
        xg::XG* index = new xg::XG();
        
        // Load it
        index->load(input);
        
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

