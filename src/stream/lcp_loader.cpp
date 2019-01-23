/**
 * \file gcsa_loader.cpp
 * Defines loaders for a GCSA LCPArray from stream files.
 */

#include "registry.hpp"

#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>

namespace vg {

namespace stream {

using namespace std;

/// Explain how to load LCP data from a series of messages viewed as a stream
static Loader<gcsa::LCPArray> lcp_loader("LCP", wrap_stream_loader([](istream& input) -> void* {
    // Allocate an LCPArray
    gcsa::LCPArray* index = new gcsa::LCPArray();
    
    // Load it
    index->load(input);
    
    // Return it so the caller owns it.
    return (void*) index;
}), wrap_stream_saver([](const void* index_void, ostream& output) {
    // Cast to LCP and serialize to the stream.
    sdsl::serialize(*(const gcsa::LCPArray*) index_void, output);
}));

}

}

