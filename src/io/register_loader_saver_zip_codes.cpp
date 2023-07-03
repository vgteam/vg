/**
 * \file register_loader_saver_zip_codes.cpp
 * Defines IO for an ZipCode index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_zip_codes.hpp"

#include "../zip_code.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_zip_codes() {

    Registry::register_bare_loader_saver_with_magic_and_filename<ZipCodeCollection>("ZIPCODES", ZipCodeCollection::get_magic_number_as_string(), 
    [](istream& input, const string& filename) -> void* {
        // Allocate an index and hand it the stream
        ZipCodeCollection* zipcodes = new ZipCodeCollection();
        if (!filename.empty()) {
            ifstream in (filename);
            zipcodes->deserialize(in);
        } else {
            zipcodes->deserialize(input);
        }
        
        // Return it so the caller owns it.
        return (void*) zipcodes;
    }, 
    [](const void* index_void, ostream& output) {
        // Cast to SnarlDistanceIndex and serialize to the stream.
        assert(index_void != nullptr);
        static_cast<const ZipCodeCollection*>(index_void)->serialize(output);
    });
}

}

}

