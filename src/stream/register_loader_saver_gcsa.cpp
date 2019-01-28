/**
 * \file register_loader_saver_gcsa.cpp
 * Defines loaders for a GCSA index from stream files.
 */

#include "registry.hpp"
#include "register_loader_saver_gcsa.hpp"

#include <gcsa/gcsa.h>

namespace vg {

namespace stream {

using namespace std;

void register_loader_saver_gcsa() {
    Registry::register_loader_saver<gcsa::GCSA>("GCSA", wrap_stream_loader([](istream& input) -> void* {
        // Allocate a GCSA
        gcsa::GCSA* index = new gcsa::GCSA();
        
        // Load it
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }), wrap_stream_saver([](const void* index_void, ostream& output) {
        // Cast to GCSA and serialize to the stream.
        assert(index_void != nullptr);
        ((const gcsa::GCSA*) index_void)->serialize(output);
    }));
}

}

}

