/**
 * \file register_loader_saver_gcsa.cpp
 * Defines IO for a GCSA index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gcsa.hpp"

#include <gcsa/gcsa.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gcsa() {
    std::uint32_t magic_number = gcsa::GCSAHeader::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gcsa::GCSA>("GCSA", magic_string, [](istream& input) -> void* {
        // Allocate a GCSA
        gcsa::GCSA* index = new gcsa::GCSA();
        
        // Load it. In case of a failure, this will:
        // * Print an error message if sanity checks fail.
        // * Fail silently if reading the input fails.
        // The exceptions are derived from std::runtime_error.
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to GCSA and serialize to the stream.
        assert(index_void != nullptr);
        ((const gcsa::GCSA*) index_void)->serialize(output);
    });
}

}

}

