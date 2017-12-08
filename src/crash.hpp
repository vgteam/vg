#ifndef VG_CRASH_HPP_INCLUDED
#define VG_CRASH_HPP_INCLUDED

#ifndef _POSIX_C_SOURCE
    #define _POSIX_C_SOURCE 199309L
    #include <signal.h>
    #undef _POSIX_C_SOURCE
#else
    #include <signal.h>
#endif

namespace vg {

using namespace std;

// Emit a stack trace when something bad happens. Add as a signal handler with sigaction.
void emit_stacktrace(int signalNumber, siginfo_t *signalInfo, void *signalContext);

}
#endif
