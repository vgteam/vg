/**
 * \file register_loader_params_json.cpp
 * Defines IO for a VG graph from stream files of Graph objects.
 */

#include <vg/io/registry.hpp>
#include "register_loader_params_json.hpp"


namespace vg {

namespace io {

using namespace std;
using namespace vg::io;

void register_loader_params_json() {
    Registry::register_loader<std::string>("PARAMS_JSON", wrap_bare_loader([](const std::istream& stream) -> void* {
        // Read the whole stream with an iterator. See <https://stackoverflow.com/a/3203502>.
        return new std::string(std::istreambuf_iterator<char>(stream), {});
    });
}

}

}

