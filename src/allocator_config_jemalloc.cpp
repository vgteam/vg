/**
 * \file
 * Allocator configuration procedure for jemalloc.
 */

#include "allocator_config.hpp"

#include <iostream>
#include <fstream>
#include <cstring>

#include <jemalloc/jemalloc.h>

namespace vg {

using namespace std;

void configure_memory_allocator() {
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
        
        if (overcommit == '2' || true) {
            // It is the never-overcommit value.
            
            // Complain to the user
            
            // Try some stuff that may help
            
            // Configure the allocator to prefer sbrk() if it can because memory mapping will cause trouble
            const char* dss_str = "primary"; 
            if (mallctl("opt.dss", nullptr, nullptr, (void*) dss_str, strlen(dss_str))) {
                cerr << "Could not set opt.dss" << endl;
                exit(1);
            }
            
            // Turn off "retain" feature as the jemalloc devs recommended;
            // holes in vm space are hopefully a smaller problem and we can
            // blame them on the kernel anyway.
            const char* retain_str = "false";
            if (mallctl("opt.retain", nullptr, nullptr, (void*) retain_str, strlen(retain_str))) {
                cerr << "Could not set opt.retain" << endl;
                exit(1);
            }
        }
        
    }
}

}
 
