/**
 * \file register_loader_saver_gbz.cpp
 * Defines IO for GBZ from stream files.
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gbz.hpp"

#include <gbwtgraph/gbz.h>

namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_saver_gbz() {
    std::uint32_t magic_number = gbwtgraph::GBZ::Header::TAG;
    std::string magic_string(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));

    Registry::register_bare_loader_saver_with_magic<gbwtgraph::GBZ>("GBZ", magic_string, [](std::istream& input) -> void* {
        gbwtgraph::GBZ* result = new gbwtgraph::GBZ();
        // Load it. In case of a failure, this will:
        // * Throw an exception if sanity checks fail.
        // * Throw an exception or fail silently if reading the input fails.
        // The exceptions are derived from std::runtime_error.
        result->simple_sds_load(input);
        return reinterpret_cast<void*>(result);
    }, [](const void* gbz_void, std::ostream& output) {
        assert(gbz_void != nullptr);
        const gbwtgraph::GBZ* gbz = reinterpret_cast<const gbwtgraph::GBZ*>(gbz_void);
        gbz->simple_sds_serialize(output);
    });
}

}

}
