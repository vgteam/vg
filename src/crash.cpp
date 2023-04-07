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

// We need realpath stuff
#include <limits.h>

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

// Where vg issues should be reported
const char* ISSUE_URL = "https://github.com/vgteam/vg/issues/new/choose";

/// Make a horizontal line for delimiting error info.
static void draw_br() {
    for (size_t i = 0; i < 20; i++) {
        std::cerr << "━";
    }
    std::cerr << std::endl;
}

/// Print an OSC-8 link start sequence to standard error, if it is a terminal,
/// pointing to the given vg source file path.
static void start_vg_link(const std::string& file_path, int line) {
    if (!isatty(fileno(stderr))) {
        return;
    }
    
    std::string url_protocol = "file";
    std::string url_host;
    std::string url_path;
    
    char real_path_buffer[PATH_MAX + 1];
    char* abspath = realpath(file_path.c_str(), real_path_buffer);
    if (abspath != nullptr) {
        // File exists to link to!
        url_path = abspath;
        
        // The link probably needs a hostname
        char host_buffer[HOST_NAME_MAX + 1];
        if (gethostname(host_buffer, HOST_NAME_MAX) == 0) {
            url_host = host_buffer;
        }
        // And we have to pick a protocol depending on if we are local or not
        if (getenv("SSH_TTY") != nullptr) {
            // We are probably developing over SSH, so link files over SFTP.
            url_protocol = "sftp";
        }
    } else {
        // File doesn't exist relative to here. Link to Github.
        url_path = "/vgteam/vg/blob/" + Version::get_version() + "/" + file_path + "#L" + std::to_string(line);
        url_protocol = "https";
        url_host = "github.com";
    }
    std::cerr << "\e]8;;" << url_protocol << "://" << url_host << url_path << "\e\\";
}

/// Print an OSC-8 link start sequence to standard error, if it is a terminal,
/// pointing to the given URL.
static void start_link(const std::string& url) {
    if (!isatty(fileno(stderr))) {
        return;
    }
    
    std::cerr << "\e]8;;" << url << "\e\\";
}

/// Print an OSC-8 link end sequence to standard error, if it is a terminal.
static void stop_link() {
    if (!isatty(fileno(stderr))) {
        return;
    }
    std::cerr << "\e]8;;\e\\";
}

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
        // Determine where to save trace files
        std::string temp;
        char* tmpdir = getenv("TMPDIR");
        if (tmpdir) {
            temp += tmpdir;
        } else {
            temp += "/tmp";
        }
        temp += "/vg_crash_XXXXXX";
        char* tempDir = mkdtemp((char*)temp.c_str());
        dirName = tempDir;
        tempStream.open(dirName+ "/stacktrace.txt");
        out = &tempStream;
    }
    
    draw_br();
    
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
        #if (defined(__arm64__) || defined(__aarch64__))
            *out << "Stack traces are not supported on ARM Macs yet" << endl;
        #else
            // macOS does it this way on x86-64
            ip = (void*)context->uc_mcontext->__ss.__rip;
            bp = (void**)context->uc_mcontext->__ss.__rbp;
            *out << "Caught signal " << signalNumber << " raised at address " << ip << endl;
            // Do our own tracing because backtrace doesn't really work on all platforms.
            stacktrace_manually(*out, signalNumber, ip, bp);
        #endif
    #elif __x86_64__
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

    // Use OSC-8 to link the user to their destination.
    cerr << "ERROR: Signal "<< signalNumber << " occurred. VG has crashed. ";
    start_link(ISSUE_URL);
    cerr << "Visit ";
    cerr << ISSUE_URL;
    cerr << " to report a bug.";
    stop_link();
    cerr << endl;
    if (fullTrace) {
        cerr << "Please include this entire error log in your bug report!" << endl; 
    } else {
        // Print path for stack trace file
        cerr << "Stack trace path: "<< dirName << "/stacktrace.txt" << endl;
        cerr << "Please include the stack trace file in your bug report!" << endl;
    }
    draw_br();
    // Make sure to exit with the right code
    exit(signalNumber + 128);
}

void enable_crash_handling() {
    // Set up stack trace support
    if (getenv(var) != nullptr) {
        if (strcmp(getenv(var), "0") == 0) {
            // if VG_FULL_TRACEBACK env var is set to 0
            fullTrace = false;
        } else {
            // if VG_FULL_TRACEBACK env var is set to anything else
            fullTrace = true;
        }
    }
    else {
        // if VG_FULL_TRACEBACK env var is not set
        fullTrace = true;
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
    sigaction(SIGFPE, &sig_config, nullptr);
    
    // We don't set_terminate for aborts because we still want the standard
    // library's message about what the exception was.
}

thread_local std::string stored_crash_context;

void set_crash_context(const std::string& message) {
    stored_crash_context = message;
}

void clear_crash_context() {
    stored_crash_context.clear();
}

void with_exception_handling(const std::function<void(void)>& body) {
    try {
        body();
    } catch(const std::exception& ex) {
        report_exception(ex); 
    }
}

void report_exception(const std::exception& ex) {
    std::cerr << "Unhandled exception: " << ex.what() << std::endl;
    if (!stored_crash_context.empty()) {
        std::cerr << "Exception context: " << stored_crash_context << std::endl;
    }
    abort();
}

void crash_unless_impl(bool condition, const std::string& condition_string, const std::string& file, int line, const std::string& function) {
    if (condition) {
        // Nothing is wrong!
        return;
    }
    std::cerr << std::endl << std::endl;
    draw_br();
    std::cerr << "VG has crashed because " << condition_string << " is false." << std::endl;
    std::cerr << "Problem is at ";
    // Use OSC-8 to link our files if we can.
    start_vg_link(file, line);
    std::cerr << file;
    stop_link();
    std::cerr << ":" << line << " in " << function << "." << std::endl;
    if (!stored_crash_context.empty()) {
        std::cerr << "This is in the context of: " << stored_crash_context << std::endl;
    }
    abort();
}

}
