// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

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

    Registry::register_bare_loader_saver<gbwt::FastLocate>("R-INDEX", [](istream& input) -> void* {
        // Allocate an r-index
        gbwt::FastLocate* index = new gbwt::FastLocate();

        // Load it
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

