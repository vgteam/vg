/**
 * \file register_loader_saver_xg.cpp
 * Defines IO for an XG index from stream files.
 */

#include "registry.hpp"
#include "register_loader_saver_gcsa.hpp"

#include "../xg.hpp"

namespace vg {

namespace stream {

using namespace std;

void register_loader_saver_xg() {
    Registry::register_bare_loader_saver<xg::XG>("XG", [](istream& input) -> void* {
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

