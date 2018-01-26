#include "crash.hpp"

// Pull in backward-cpp
#ifndef __APPLE__
    #define BACKWARD_HAS_DW 1
#endif
#include <backward.hpp>


namespace vg {

void enable_crash_handling() {
    // Set up stack trace support

    // Use backward-cpp
    static backward::SignalHandling signal_handling_manager;

    // We don't set_terminate for aborts because we still want the standard
    // library's message about what the exception was.
}

}
