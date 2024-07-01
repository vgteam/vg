#ifndef VG_ALLOCATOR_CONFIG_HPP_INCLUDED
#define VG_ALLOCATOR_CONFIG_HPP_INCLUDED

/**
 * \file
 * Allocator configuration header. Used with either
 * allocator_config_jemalloc.cpp or allocator_config_system.cpp as appropriate
 * for the build.
 *
 * Contains startup functions and functions to manipulate memory profiling, if available.
 */
 
namespace vg {

/**
 * Interface for working with the memory allocator that is compiled into the build.
 */
struct AllocatorConfig {

    /**
     * If using a non-system memory allocator, initialize it to a safe
     * configuration in this runtime environment.
     */
    static void configure();

    /**
     * Turn memory profiling on or off, if available in the allocator.
     */
    static void set_profiling(bool should_profile);

    /**
     * Dump a memory profiling snapshot, if available in the allocator.
     */
    static void snapshot();

};

}
 
#endif
