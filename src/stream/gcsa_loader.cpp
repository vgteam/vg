/**
 * \file gcsa_loader.cpp
 * Defines loaders for a GCSA index from stream files.
 */

#include "registry.hpp"

#include <gcsa/gcsa.h>

namespace vg {

namespace stream {

using namespace std;

/// Explain how to load GCSA data from a series of messages viewed as a stream
static Loader<gcsa::GCSA> gcsa_loader("GCSA", wrap_stream_loader([](istream& input) -> void* {
    // Allocate a GCSA
    gcsa::GCSA* index = new gcsa::GCSA();
    
    // Load it
    index->load(input);
    
    // Return it so the caller owns it.
    return (void*) index;
}), wrap_stream_saver([](const void* index_void, ostream& output) {
    // Cast to GCSA and serialize to the stream.
    sdsl::serialize(*(const gcsa::GCSA*) index_void, output);
}));

}

}

