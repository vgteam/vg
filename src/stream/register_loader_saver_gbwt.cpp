/**
 * \file register_loader_saver_gbwt.cpp
 * Defines IO for a GBWT index from stream files.
 */

#include "registry.hpp"
#include "register_loader_saver_gcsa.hpp"

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>

namespace vg {

namespace stream {

using namespace std;

void register_loader_saver_gbwt() {
    // GBWT and DynamicGBWT can both load/save the same format.
    // TODO: Should DynamicGBWT be its own file here?

    Registry::register_bare_loader_saver<gbwt::GBWT>("GBWT", [](istream& input) -> void* {
        // Allocate a GBWT
        gbwt::GBWT* index = new gbwt::GBWT();
        
        // Load it
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to GBWT and serialize to the stream.
        assert(index_void != nullptr);
        ((const gbwt::GBWT*) index_void)->serialize(output);
    });
    
    Registry::register_bare_loader_saver<gbwt::DynamicGBWT>("GBWT", [](istream& input) -> void* {
        // Allocate a DynamicGBWT
        gbwt::DynamicGBWT* index = new gbwt::DynamicGBWT();
        
        // Load it
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to DynamicGBWT and serialize to the stream.
        assert(index_void != nullptr);
        ((const gbwt::DynamicGBWT*) index_void)->serialize(output);
    });
}

}

}

