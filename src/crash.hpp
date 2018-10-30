#ifndef VG_CRASH_HPP_INCLUDED
#define VG_CRASH_HPP_INCLUDED

/**
 * \file crash.hpp
 *
 * Implementation for crash handling to create a stack trace when VG crashes.
 * To use the crash handling system, call enable_crash_handling() early on in the program.
 * When a crash occurs, you will recieve an error message with the path to the 
 * stack trace file. To get the full stack trace on standard error, you need to
 * set the environment variable 'VG_FULL_TRACEBACK=1'. 
 *
 */

namespace vg {

/// Main should call this to turn on our stack tracing support.
void enable_crash_handling();

}
#endif
