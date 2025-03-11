/**
 * \file register_loader_saver_r_index.cpp
 * Defines IO for an r-index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_r_index.hpp"

#include <gbwt/fast_locate.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_r_index() {
    std::uint32_t magic_number = gbwt::FastLocate::Header::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gbwt::FastLocate>("R-INDEX", magic_string, [](istream& input) -> void* {
        // Allocate an r-index
        gbwt::FastLocate* index = new gbwt::FastLocate();

        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Fail silently if reading the input fails.
        // The exceptions are derived from std::runtime_error.
        index->load(input);

        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to r-index and serialize to the stream.
        assert(index_void != nullptr);
        ((const gbwt::FastLocate*) index_void)->serialize(output);
    });
}

}

}

