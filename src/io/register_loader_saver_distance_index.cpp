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

    bdsg::SnarlDistanceIndex empty;
    Registry::register_bare_loader_saver_with_magic_and_filename<SnarlDistanceIndex>("DISTANCE2", empty.get_prefix(), 
    [](istream& input, const string& filename) -> void* {
        // Allocate an index and hand it the stream
        SnarlDistanceIndex* index = new SnarlDistanceIndex();
        if (!filename.empty()) {
            index->deserialize(filename);
        } else {
            index->deserialize(input);
        }
        
        // Return it so the caller owns it.
        return (void*) index;
    }, 
    [](const void* index_void, ostream& output) {
        // Cast to SnarlDistanceIndex and serialize to the stream.
        assert(index_void != nullptr);
        throw std::runtime_error( "warning [vpkg::save<DistanceIndex>]: save the distance index directly with serialize() instead of with vpkg");

        //((const SnarlDistanceIndex*) index_void)->serialize(output);
    });
}

}

}

