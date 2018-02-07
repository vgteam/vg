#ifndef VG_CRASH_HPP_INCLUDED
#define VG_CRASH_HPP_INCLUDED

// You will need to #define _POSIX_C_SOURCE 199309L or greater for this header
// to work, because the types it uses from signal.h are behind that macro.
#include <signal.h>

namespace vg {

using namespace std;

/// Main should call this to turn on our stack tracing support.
void enable_crash_handling();

}
#endif
