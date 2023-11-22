/**
 * \file
 * Allocator configuration procedure for jemalloc.
 */

#include "allocator_config.hpp"

#include <iostream>
#include <fstream>
#include <cstring>

#include <jemalloc/jemalloc.h>

extern "C" {
    // Hackily define symbols that jemalloc actually exports.
    // Somehow it gets a "je_" prefix on these relative to what's in its
    // source.
    // They're also all "local" symbols in the dynamic jemalloc library,
    // meaning we can't link them from outside the library; we need to use
    // static jemalloc if we intend to access these from here.
    
    // We use int here but really this takes an enum type.
    bool je_extent_dss_prec_set(int dss_prec);
    
    // This is really the last enum value
    int dss_prec_limit = 3;

    // These are the globals used to store the human-readable dss priority in
    // addition to what the function controls.
    extern const char *je_opt_dss;
    extern const char *je_dss_prec_names[];
    
    extern bool je_opt_retain;
}

// Stringifier we need for jemalloc from its docs
#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

namespace vg {

using namespace std;

void AllocatorConfig::configure() {
    // TODO: this is going to allocate when we don't really maybe want to. But
    // the dynamic linker also allocated; we have to hope we don't upset any
    // existing jemalloc stuff.
    ifstream procfile("/proc/sys/vm/overcommit_memory");
    if (procfile) {
        // We're actually on a Linux system with an overcommit setting.
        // TODO: Can it be changed on Mac?
        
        // We need to work around jemalloc's propensity to run out of memory
        // mappings and fail to allocate, when overcommit is disabled and the
        // number of distinct mappings is capped. See <https://github.com/jemalloc/jemalloc/issues/1328>
        
        // Read the setting
        char overcommit;
        procfile >> overcommit;
        
        if (overcommit == '2') {
            // It is the never-overcommit value.
            
            // Complain to the user
            cerr << "vg [warning]: System's vm.overcommit_memory setting is 2 (never overcommit). "
                << "vg does not work well under these conditions; you may appear to run out of memory with plenty of memory left. "
                << "Attempting to unsafely reconfigure jemalloc to deal better with this situation." << endl;
            
            // Try some stuff that may help
            
            // Configure the allocator to prefer sbrk() if it can because memory mapping will cause trouble
            const char* dss_str = "primary";
            size_t dss_str_len = strlen(dss_str);
            
            bool match = false;
            // Redo the dss_prec loop from jemalloc: <https://github.com/jemalloc/jemalloc/blob/83f3294027952710f35014cff1cffd51f281d785/src/jemalloc.c#L1194-L1208>
            // This should cover newly created arenas.
            for (int i = 0; i < dss_prec_limit; i++) {
                if (strncmp(je_dss_prec_names[i], dss_str, dss_str_len) == 0) {
                    if (je_extent_dss_prec_set(i)) {
                        cerr << "Could not reconfigure jemalloc dss_prec" << endl;
                        exit(1);
                    } else {
                        je_opt_dss = je_dss_prec_names[i];
                        match = true;
                        break;
                    }
                }
            }
            if (!match) {
                cerr << "Could not find jemalloc dss_prec of " << dss_str << endl;
                exit(1);
            }
            // Then fix up all existing arenas (without allocating?)
            // To write these string parameters we need to copy a pointer into place, not a value
            const char** dss_str_location = &dss_str; 
            auto mallctl_result = mallctl("arena." STRINGIFY(MALLCTL_ARENAS_ALL) ".dss", nullptr, nullptr, (void*) dss_str_location, sizeof(dss_str_location));
            if (mallctl_result) {
                cerr << "Could not set dss priority on existing jemalloc arenas: " << strerror(mallctl_result) << endl;
                exit(1);
            }
            
            // Finally, make the opt_retain flag be off.
            // This seems most likely to upset jemalloc because it changes the semantics of some of its internal fields.
            je_opt_retain = false;
        }
        
    }
}

void AllocatorConfig::set_profiling(bool should_profile) {
    // Send the bool right into jemalloc's profiling-is-active flag.
    //
    // You need to start vg with something like
    // MALLOC_CONF="prof_active:false,prof:true" for this to be useful.
    auto mallctl_result = mallctl("prof.active", nullptr, nullptr, &should_profile, sizeof(should_profile));
    if (mallctl_result) {
        std::cerr << "Could not set profiling to " << should_profile << ": " << strerror(mallctl_result) << std::endl;
        exit(1);
    }
}

void AllocatorConfig::snapshot() {
    // Ask to dump a profile now.
    //
    // You need to start vg with something like
    // MALLOC_CONF="prof_prefix:jeprof.out" for this to have a filename to go
    // to.
    auto mallctl_result = mallctl("prof.dump", NULL, NULL, NULL, 0);
    if (mallctl_result) {
        std::cerr << "Could not dump profile: " << strerror(mallctl_result) << std::endl;
        exit(1);
    }
}

}
 
