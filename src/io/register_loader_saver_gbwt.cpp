/**
 * \file register_loader_saver_gbwt.cpp
 * Defines IO for a GBWT index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbwt.hpp"

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbwt() {
    // GBWT and DynamicGBWT can both load/save the same format.
    std::uint32_t magic_number = gbwt::GBWTHeader::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gbwt::GBWT>("GBWT", magic_string, [](istream& input) -> void* {
        // Allocate a GBWT
        gbwt::GBWT* index = new gbwt::GBWT();

        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Throw an exception or fail silently if reading a Simple-SDS input fails.
        // * Exit with std::exit() or fail silently if reading an SDSL input fails.
        // The exceptions are derived from std::runtime_error.
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to GBWT and serialize to the stream.
        assert(index_void != nullptr);
        ((const gbwt::GBWT*) index_void)->simple_sds_serialize(output);
    });
    
    Registry::register_bare_loader_saver_with_magic<gbwt::DynamicGBWT>("GBWT", magic_string, [](istream& input) -> void* {
        // Allocate a DynamicGBWT
        gbwt::DynamicGBWT* index = new gbwt::DynamicGBWT();

        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Throw an exception or fail silently if reading a Simple-SDS input fails.
        // * Exit with std::exit() or fail silently if reading an SDSL input fails.
        // The exceptions are derived from std::runtime_error.
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to DynamicGBWT and serialize to the stream.
        assert(index_void != nullptr);
        ((const gbwt::DynamicGBWT*) index_void)->simple_sds_serialize(output);
    });
}

}

}

