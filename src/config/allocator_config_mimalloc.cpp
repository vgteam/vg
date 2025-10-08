/**
 * \file
 * Allocator configuration procedure for mimalloc.
 */

#include "allocator_config.hpp"

#include <mimalloc.h>

namespace vg {

using namespace std;

void AllocatorConfig::configure() {
    // mimalloc's default configuration is fine, so do nothing.
}

void AllocatorConfig::set_profiling(bool should_profile) {
    // mimalloc doesn't have builtin profiling, so do nothing.
}

void AllocatorConfig::snapshot() {
    // mimalloc doesn't have builtin profiling, so do nothing.
}

}
 
