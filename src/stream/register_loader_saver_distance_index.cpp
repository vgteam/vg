/**
 * \file register_loader_saver_distance_index.cpp
 * Defines IO for an XG index from stream files.
 */

#include "registry.hpp"
#include "register_loader_saver_distance_index.hpp"

#include "../distance.hpp"

namespace vg {

namespace stream {

using namespace std;

void register_loader_saver_distance_index() {
    Registry::register_bare_loader_saver<DistanceIndex>("DISTANCE", [](istream& input) -> void* {
        // Allocate an index and hand it the stream
        DistanceIndex* index = new DistanceIndex(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to DistanceIndex and serialize to the stream.
        assert(index_void != nullptr);
        ((const DistanceIndex*) index_void)->serialize(output);
    });
}

}

}

