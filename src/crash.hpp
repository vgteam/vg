#ifndef VG_CRASH_HPP_INCLUDED
#define VG_CRASH_HPP_INCLUDED

/**
 * \file crash.hpp
 *
 * Implementation for crash handling to create a stack trace when VG crashes.
 * To use the crash handling system, call enable_crash_handling() early on in the program.
 * When a crash occurs, you will recieve an error message with the stack trace.
 * To get just a filename, you need to set the environment variable
 * 'VG_FULL_TRACEBACK=0'. 
 *
 */

#include <functional>
#include <string>

namespace vg {

/// Main should call this to turn on our stack tracing support.
void enable_crash_handling();

/// User code should call this when it has context for a failure in its thread.
void set_crash_context(const std::string& message);

/// User code should call this when it wants to clear context for a failure in its thread.
void clear_crash_context();

/// User code should call this to get all its exceptions handled.
void with_exception_handling(const std::function<void(void)>& body);

/// User code should call this if it catches an exception it doesn't know what
/// to do with.
void report_exception(const std::exception& ex);

/// User code should call this instead of assert
#define crash_unless(condition) crash_unless_impl((condition), #condition, __FILE__, __LINE__, __func__); 

/// crash_unless calls into this function for a real implementation.
void crash_unless_impl(bool condition, const std::string& condition_string, const std::string& file, int line, const std::string& function);


}
#endif
