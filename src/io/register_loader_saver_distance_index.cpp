/**
 * \file register_loader_saver_distance_index.cpp
 * Defines IO for an XG index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_distance_index.hpp"

#include "../snarl_distance_index.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_distance_index() {
    // The distance index header is just a text string. We need to make sure
    // this looks like a bare distance index file if we are going to load it
    // without type-tagged message deserialization.
    Registry::register_bare_loader_saver_with_magic<SnarlDistanceIndex>("DISTANCE2", "distance index", [](istream& input) -> void* {
        // Allocate an index and hand it the stream
        SnarlDistanceIndex* index;
        index->deserialize_members(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to SnarlDistanceIndex and serialize to the stream.
        assert(index_void != nullptr);
        ((const SnarlDistanceIndex*) index_void)->serialize_members(output);
    });
}

}

}

