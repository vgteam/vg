// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

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
    Registry::register_bare_loader_saver<gcsa::GCSA>("GCSA", [](istream& input) -> void* {
        // Allocate a GCSA
        gcsa::GCSA* index = new gcsa::GCSA();
        
        // Load it
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

