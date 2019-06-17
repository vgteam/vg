/**
 * \file register_loader_saver_distance_index.cpp
 * Defines IO for an xg::XG index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_distance_index.hpp"

#include "../min_distance.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_distance_index() {
    Registry::register_bare_loader_saver<MinimumDistanceIndex>("DISTANCE", [](istream& input) -> void* {
        // Allocate an index and hand it the stream
        MinimumDistanceIndex* index = new MinimumDistanceIndex(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to MinimumDistanceIndex and serialize to the stream.
        assert(index_void != nullptr);
        ((const MinimumDistanceIndex*) index_void)->serialize(output);
    });
}

}

}

