#ifdef __APPLE__
    // We can't set _POSIX_C_SOURCE on OS X because it hides the non-POSIX dladdr() with no way to get it back.
    // But we do need level 200809L features generally (for example, for EOWNERDEAD and ENOTRECOVERABLE on OS X).
    // So we do it the Apple Way.
    #define __DARWIN_C_LEVEL __DARWIN_C_FULL
#else
    // Needed for siginfo_t to be defined.
    #define _POSIX_C_SOURCE 200809L
#endif

// We need this for ucontext.h to give us useful context stuff...
#define _XOPEN_SOURCE 1

#include "crash.hpp"
#include "version.hpp"

// iostream wants this on Travis on Mac
#include <pthread.h>
#include <errno.h>

// Needed for automatic name demangling, but not all that portable
#include <cxxabi.h>
#include <execinfo.h>

// Needed to hack contexts for signal traces
 #include <ucontext.h>
 
// Needed to generate stacktraces ourselves
#include <dlfcn.h>

// We need strcmp
#include <cstring>

#include <cstdlib>
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <signal.h>

#ifndef __APPLE__
    // Pull in backward-cpp and use libdw from elfutils.
    // In theory backward-cpp can build and even backtrace on mac
    // In practice the mac port doesn't work on my machine and breaks the build on Travis.
    #define BACKWARD_HAS_DW 1
    #include <backward.hpp>
#endif

namespace vg {

// env var for getting full stack trace on cerr instead of a file path
const char* var = "VG_FULL_TRACEBACK";
// fullTrace = true means env var was set
// fullTrace = false means env var was not set
bool fullTrace = false;

void stacktrace_manually(ostream& out, int signalNumber, void* ip, void** bp) {
    // Now we compute our own stack trace, because backtrace() isn't so good on OS X.
    // We operate on the same principles as <https://stackoverflow.com/a/5426269>

    // Unfortunately this *only* seems to work even a little well on OS X; on Linux we get ip
    // wandering off relatively quickly.
    
    // Allocate a place to keep the dynamic library info for the address the stack is executing at
    Dl_info address_library;
    
    out << endl;
    out << "Next ip: " << ip << " Next bp: " << bp << endl;
    while (true) {
        // Work out where the ip is.
        if (!dladdr(ip, &address_library)) {
            // This address doesn't belong to anything!
            out << "Stack leaves code at ip=" << ip << endl;
            break;
        }

        if (address_library.dli_sname != nullptr) {
            // Make a place for the demangling function to save its status
            int status;
            
            // Do the demangling
            char* demangledName = abi::__cxa_demangle(address_library.dli_sname, NULL, NULL, &status);
            
            if(status == 0) {
                // Successfully demangled
                out << "Address " << ip << " in demangled symbol " << demangledName
                    << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_saddr))
                    << ", in library " << address_library.dli_fname
                    << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_fbase)) << endl;
                free(demangledName);
            } else {
                // Leave mangled
                out << "Address " << ip << " in mangled symbol " << address_library.dli_sname
                    << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_saddr))
                    << ", in library " << address_library.dli_fname
                    << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_fbase)) << endl;
            }
            
            #ifdef __APPLE__
                #ifdef VG_DO_ATOS
                    // Try running atos to print a line number. This can be slow so we don't do it by default.
                    stringstream command;
                    
                    command << "atos -o " << address_library.dli_fname << " -l " << address_library.dli_fbase << " " << ip;
                    out << "Running " << command.str() << "..." << endl;
                    system(command.str().c_str());
                #endif
            #endif
            
        } else {
            out << "Address " << ip << " out of symbol in library " << address_library.dli_fname << endl;
        }

        if(address_library.dli_sname != nullptr && !strcmp(address_library.dli_sname, "main")) {
            out << "Stack hit main" << endl;
            break;
        }

        if (bp != nullptr) {
            // Simulate a return
            ip = bp[1];
            bp = (void**) bp[0];
        } else {
            break;
        }
        
        out << endl;
        out << "Next ip: " << ip << " Next bp: " << bp << endl;
    }
    
    out << "Stack trace complete" << endl;
    out << endl;
    
}

/// Emit a stack trace when something bad happens. Add as a signal handler with sigaction.
void emit_stacktrace(int signalNumber, siginfo_t *signalInfo, void *signalContext) {
    ofstream tempStream;
    string dirName;
    ostream* out;

    if (fullTrace == true) {
        out = &cerr;
    } else {
        char temp[] = "/tmp/vg_crash_XXXXXX";
        char* tempDir = mkdtemp(temp);
        dirName = tempDir;
        tempStream.open(dirName+ "/stacktrace.txt");
        out = &tempStream;
    }
    
    *out << "Crash report for vg " << Version::get_short() << endl;
   
    // This holds the context that the signal came from, including registers and stuff
    ucontext_t* context = (ucontext_t*) signalContext;
   
    // TODO: This assumes x86_64
    // Fetch out the registers
    // We model IP as a pointer to void (i.e. into code)
    void* ip;
    // We model BP as an array of two things: previous BP, and previous IP.
    void** bp;
   
    #ifdef __APPLE__
        // OS X 64 bit does it this way
        ip = (void*)context->uc_mcontext->__ss.__rip;
        bp = (void**)context->uc_mcontext->__ss.__rbp;
        *out << "Caught signal " << signalNumber << " raised at address " << ip << endl;
        // Do our own tracing because backtrace doesn't really work on all platforms.
        stacktrace_manually(*out, signalNumber, ip, bp);
    #else
        // Linux 64 bit does it this way
        ip = (void*)context->uc_mcontext.gregs[REG_RIP];
        bp = (void**)context->uc_mcontext.gregs[REG_RBP];
        
        static backward::StackTrace stack_trace;
        stack_trace.load_from(ip, 32);
        static backward::Printer p;
        p.color_mode = backward::ColorMode::automatic;
        p.address = true;
        p.object = true;
        p.print(stack_trace, *out);
        tempStream.close();
    #endif

    if (fullTrace == false) {
        cerr << "ERROR: Signal "<< signalNumber << " occurred. VG has crashed. Run 'vg bugs --new' to report a bug." << endl;
        // Print path for stack trace file
        cerr << "Stack trace path: "<< dirName << "/stacktrace.txt" << endl;
    }
    // Make sure to exit with the right code
    exit(signalNumber + 128);
}

void enable_crash_handling() {
    // Set up stack trace support
    if (getenv(var) != nullptr) {
        if (strcmp(getenv(var), "1") == 0) {
            // if VG_FULL_TRACEBACK env var is set
            fullTrace = true;
        }
    }
    else {
        // if VG_FULL_TRACEBACK env var is not set
        fullTrace = false;
    }
    // backtrace() doesn't work in our Mac builds, and backward-cpp uses backtrace().
    // Do this the old-fashioned way.
    
    // We do it the cleverer sigaction way to try and make OS X backtrace not just tell us that the signal handler is being called.
    struct sigaction sig_config;
    sig_config.sa_flags = SA_SIGINFO; // Use the new API and not the old signal() API for the handler.
    sig_config.sa_sigaction = emit_stacktrace;
    sigemptyset(&sig_config.sa_mask);
 
    sigaction(SIGABRT, &sig_config, nullptr);
    sigaction(SIGSEGV, &sig_config, nullptr);
    sigaction(SIGBUS, &sig_config, nullptr);
    sigaction(SIGILL, &sig_config, nullptr);
    
    // We don't set_terminate for aborts because we still want the standard
    // library's message about what the exception was.
}

}
