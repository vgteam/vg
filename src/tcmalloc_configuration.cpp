/**
 * \file tcmalloc_configuration.cpp
 * Defines static tcmalloc configuration.
 */

#include <gperftools/malloc_extension.h>

#include <cassert>

#include "utility.hpp"

namespace vg {

using namespace std;

/**
 * A class which, when constructed, sets up the tcmalloc allocator.
 * We make a static instance in this compilation unit.
 * We link it in only when we are using tcmalloc.
 */
class TCMallocConfiguration {
public:
    TCMallocConfiguration() {
    
        // set a higher value for tcmalloc warnings
        setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);
        
        // Allow 8 MB of memory per core (i.e. per default OMP thread) for tcmalloc thread caches.
        // The default is 16 MB total, which is way too small.
        size_t tcmalloc_thread_cache_bytes = get_thread_count() * 8 * 1024 * 1024;
        // Ask tcmalloc to up our thread cache size.
        // Don't stop if this fails; it doesn't work under callgrind for example.
        MallocExtension::instance()->SetNumericProperty("tcmalloc.max_total_thread_cache_bytes", tcmalloc_thread_cache_bytes);
    }
};

/// Create a static instance that actually applies the configuration when linked in.
static TCMallocConfiguration configure_tcmalloc;

}
