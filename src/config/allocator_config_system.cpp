/**
 * \file
 * Allocator configuration procedure for the system allocator.
 */

#include "allocator_config.hpp"
#include <cstdlib>

#ifdef __GLIBC__
// We need a bunch of machinery for using glibc's malloc_info.
#include <atomic>
#include <cstdio>
#include <sstream>
#include <malloc.h>
#endif

namespace vg {

void AllocatorConfig::configure() {
    // Nothing to do! The system allocator may be slow or not, depending on the
    // system, but it isn't really configurable in any meaningful way.
}

void AllocatorConfig::set_profiling(bool should_profile) {
    // Nothing to do! There is no standard profiling interface.
}

void AllocatorConfig::snapshot() {
#ifdef __GLIBC__
    // Track snapshot number so each snapshot is distinct.
    static std::atomic<int> snapshot_number(0);
    // Make up a filename
    std::stringstream ss;
    ss << "malloc_info.";
    ss << snapshot_number.fetch_add(1);
    ss << ".xml";
    
    // Opejn the file
    FILE* dumpfile = fopen(ss.str().c_str(), "w");
    if (dumpfile) {
        // And if that worked, dump to it.
        malloc_info(0, dumpfile);
        // And close it
        fclose(dumpfile);
    }
#endif
}

}
 

