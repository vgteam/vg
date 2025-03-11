/**
 * \file register_loader_saver_gbzgraph.cpp
 * Defines IO for GBZ in a handle graph proxy wrapper from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbzgraph.hpp"

#include <gbwtgraph/gbz.h>
#include "../gbzgraph.hpp"

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbzgraph() {
    std::uint32_t magic_number = gbwtgraph::GBZ::Header::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<GBZGraph, PathHandleGraph, HandleGraph>("GBZ", magic_string, [](std::istream& input) -> void* {
        GBZGraph* result = new GBZGraph();
        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Throw an exception or fail silently if reading the input fails.
        // The exceptions are derived from std::runtime_error.
        result->gbz.simple_sds_load(input);
        return reinterpret_cast<void*>(result);
    }, [](const void* gbzgraph_void, std::ostream& output) {
        assert(gbzgraph_void != nullptr);
        const GBZGraph* gbz_graph = reinterpret_cast<const GBZGraph*>(gbzgraph_void);
        gbz_graph->gbz.simple_sds_serialize(output);
    });
}

}

}
