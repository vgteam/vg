/**
 * \file register_loader_saver_minimizer.cpp
 * Defines IO for a minimizer index from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_minimizer.hpp"

#include <gbwtgraph/minimizer.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_minimizer() {
    std::uint32_t magic_number = gbwtgraph::MinimizerHeader::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gbwtgraph::DefaultMinimizerIndex>("MinimizerIndex", magic_string, [](istream& input) -> void* {
        gbwtgraph::DefaultMinimizerIndex* index = new gbwtgraph::DefaultMinimizerIndex();
        index->deserialize(input);
        
        // Return the index so the caller owns it.
        return static_cast<void*>(index);
    }, [](const void* index_void, ostream& output) {
        assert(index_void != nullptr);
        static_cast<const gbwtgraph::DefaultMinimizerIndex*>(index_void)->serialize(output);
    });
}

}

}

