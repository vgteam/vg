#ifndef VG_CRASH_HPP_INCLUDED
#define VG_CRASH_HPP_INCLUDED

// You will need to #define _POSIX_C_SOURCE 199309L or greater for this header
// to work, because the types it uses from signal.h are behind that macro.
#include <signal.h>

namespace vg {

using namespace std;

// Emit a stack trace when something bad happens. Add as a signal handler with sigaction.
void emit_stacktrace(int signalNumber, siginfo_t *signalInfo, void *signalContext);

}
#endif
