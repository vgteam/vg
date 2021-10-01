/**
 * \file
 * Allocator configuration procedure for the system allocator.
 */

#include "allocator_config.hpp"

namespace vg {

using namespace std;

void configure_memory_allocator() {
    // Nothing to do! The system allocator may be slow or not, depending on the
    // system, but it isn't really configurable in any meaningful way.
}

}
 

