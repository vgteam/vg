/**
 * \file gcsa_loader.cpp
 * Defines loaders for GCSA types from stream files.
 */

#include "registry.hpp"

#include <gcsa/gcsa.h>

namespace vg {

namespace stream {

using namespace std;

/// Explain how to load GCSA data from a series of messages.
static Loader<gcsa::GCSA> gcsa_loader("GCSA", [&](const message_sender_function_t& for_each_message) -> void* {
    // TODO: Implement
    return nullptr;
});

}

}

