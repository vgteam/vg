#ifndef VG_ALLOCATOR_CONFIG_HPP_INCLUDED
#define VG_ALLOCATOR_CONFIG_HPP_INCLUDED

/**
 * \file
 * Allocator configuration header. Used with either
 * allocator_config_jemalloc.cpp or allocator_config_system.cpp as appropriate
 * for the build.
 */
 
namespace vg {

/**
 * If using a non-system memory allocator, initialize it to a safe configuration in this runtime environment.
 */
void configure_memory_allocator();

}
 
#endif
