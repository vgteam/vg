/**
 * \file register_loader_saver_lcp.cpp
 * Defines IO for a GCSA LCPArray from stream files.
 */

#include <vg/io/registry.hpp>

#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_lcp() {
    std::uint32_t magic_number = gcsa::LCPHeader::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gcsa::LCPArray>("LCP", magic_string, [](istream& input) -> void* {
        // Allocate an LCPArray
        gcsa::LCPArray* index = new gcsa::LCPArray();
        
        // Load it
        index->load(input);
        
        // Return it so the caller owns it.
        return (void*) index;
    }, [](const void* index_void, ostream& output) {
        // Cast to LCP and serialize to the stream.
        ((const gcsa::LCPArray*) index_void)->serialize(output);
    });

}

}

}

