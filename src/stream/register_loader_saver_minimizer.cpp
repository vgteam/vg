/**
 * \file register_loader_saver_minimizer.cpp
 * Defines IO for a minimizer index from stream files.
 */

#include "registry.hpp"
#include "register_loader_saver_minimizer.hpp"

#include "../minimizer.hpp"

namespace vg {

namespace stream {

void register_loader_saver_minimizer() {
    Registry::register_bare_loader_saver<MinimizerIndex>("MinimizerIndex", [](istream& input) -> void* {
        MinimizerIndex* index = new MinimizerIndex();
        index->load(input);
        
        // Return the index so the caller owns it.
        return static_cast<void*>(index);
    }, [](const void* index_void, ostream& output) {
        assert(index_void != nullptr);
        static_cast<const MinimizerIndex*>(index_void)->serialize(output);
    });
}

}

}

