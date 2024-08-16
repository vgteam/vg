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
#include <sstream>
#include <signal.h>
#include <typeinfo>

#if !(defined(__APPLE__) && defined(__x86_64__))
    #ifndef __APPLE__
        // Use libdw from elfutils.
        #define BACKWARD_HAS_DW 1
    #endif
    // In theory backward-cpp can build and even backtrace on mac
    #include <backward.hpp>
#endif

#include <omp.h>

namespace vg {

/// Each thread stores a string of its crash context locally for exception handling
thread_local std::string stored_crash_context;

// We also store context data statically for signal handling. This needs OMP.

/// How many chartacters of context do we store statically?
constexpr static size_t CONTEXT_BUFFER_SIZE = 256;
/// How many threads do we store static context data for?
constexpr static size_t CONTEXT_BUFFER_COUNT = 256;
/// Stores not-always-null-terminated context data. The compiler automatically
/// initializes this to nulls.
static char context_buffer[CONTEXT_BUFFER_COUNT][CONTEXT_BUFFER_SIZE];

void set_crash_context(const std::string& message) {
    // Store locally
    stored_crash_context = message;

    size_t thread_num = omp_get_thread_num();
    if (thread_num < CONTEXT_BUFFER_COUNT) {
        // Store for other threads.
        strncpy(context_buffer[thread_num], message.c_str(), CONTEXT_BUFFER_SIZE);
    }
}

void clear_crash_context() {
    // Clear locally
    stored_crash_context.clear();

    size_t thread_num = omp_get_thread_num();
    if (thread_num < CONTEXT_BUFFER_COUNT) {
        // Clear for other threads
        context_buffer[thread_num][0] = '\0';
    }
}

/**
 * Log all stored crash contexts to the given stream.
 * 
 * Will produce undefined string values if the threads in question update their
 * contexts at the same time.
 */
static void dump_crash_contexts(std::ostream& out) {
    out << "Context dump:" << std::endl;
    // We need to copy to a local buffer because the other thread may still be running!
    char local_buffer[CONTEXT_BUFFER_SIZE];
    size_t threads_with_context = 0;
    for (size_t i = 0; i < CONTEXT_BUFFER_COUNT; i++) {
        strncpy(local_buffer, context_buffer[i], CONTEXT_BUFFER_SIZE);
        if (local_buffer[0] != '\0') {
            // Somebody wrote something here and never cleared it.
            local_buffer[CONTEXT_BUFFER_SIZE - 1] = '\0';
            out << "\tThread " << i << ": " << local_buffer << std::endl;
            threads_with_context++;
        }
    }
    out << "Found " << threads_with_context << " threads with context." << std::endl;
}

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
        
        size_t host_length_limit;
        #if defined(HOST_NAME_MAX)
            host_length_limit = HOST_NAME_MAX;
        #elif defined(_POSIX_HOST_NAME_MAX)
            host_length_limit = _POSIX_HOST_NAME_MAX;
        #else
            host_length_limit = 256;
        #endif
        
        // The link probably needs a hostname
        char host_buffer[host_length_limit + 1];
        if (gethostname(host_buffer, host_length_limit) == 0) {
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

// Report a loaded library location, or an actual source file location if we can get it
// If we need to do supplemental command line command lookups of source lines that backward-cpp can't do, do those too.
// Does not include a trailing newline.
void report_library(ostream& out, Dl_info& address_library, void* ip) {
    #ifdef __APPLE__
        // Try running atos to print a line number. This can be slow so we don't do it by default.
        stringstream command;
        
        command << "atos -o " << address_library.dli_fname << " -l " << address_library.dli_fbase << " " << ip;
        
        FILE* command_pipe = popen(command.str().c_str(), "r");
        if (command_pipe != NULL) {
            // We started the command

            // Read the result. May or may not actually work, but if nothing is read it returns 0.
            char result_buffer[1024];
            size_t bytes_read = fread(result_buffer, 1, 1023, command_pipe);
            while (bytes_read != 0 && result_buffer[bytes_read - 1] == '\n') {
                // Strip off trailing newlines
                bytes_read--;
            }
            // Add null terminator.
            result_buffer[bytes_read] = 0;

            // Dump any extra bytes so we can wait on the command.
            while (fgetc(command_pipe) != EOF) {
                // Do nothing
            }

            if (pclose(command_pipe) == 0) {
                // The command ducceeded. Report what it said and the library path.
                out << result_buffer << " in " << address_library.dli_fname << " loaded at " << address_library.dli_fbase;
                return;
            }
        }
    #endif

    // If we don't quit early, just talk about the library.
    out << "Library " << address_library.dli_fname << " loaded at " << address_library.dli_fbase;
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
        } else {
            out << "Address " << ip << " out of symbol in library " << address_library.dli_fname << endl;
        }

        out << "\t";
        report_library(out, address_library, ip);
        out << std::endl;

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

    // See
    // <https://github.com/bombela/backward-cpp/blob/65a769ffe77cf9d759d801bc792ac56af8e911a3/backward.hpp#L4239>
    // for how to decode this on different platforms.
    
    #if defined(__APPLE__) && defined(__x86_64__) 
        // On x86-64 Mac we do a manual stack trace.
        // We model IP as a pointer to void, into the code(?)
        void* ip = (void*)context->uc_mcontext->__ss.__rip;
        // We model BP as an array of two things: previous BP, and previous IP.
        void** bp = (void**)context->uc_mcontext->__ss.__rbp;
        *out << "Caught signal " << signalNumber << " raised at address " << ip << endl;
        // Do our own tracing because backtrace doesn't really work on all platforms.
        stacktrace_manually(*out, signalNumber, ip, bp);
    #else
        // Everywhere else we know of, we try backward-cpp.
        // TODO: For some reason we don't need bp?
        void* ip = nullptr;

        #if defined(__APPLE__)
            // Mac (not x86_64)
            #if (defined(__arm64__) || defined(__aarch64__))
                // Arm Mac does it this way
                ip = (void*)context->uc_mcontext->__ss.__pc;
            #endif
        #else
            // Linux
            #if defined(__x86_64__)
                // Linux x86-64 does it this way
                ip = (void*)context->uc_mcontext.gregs[REG_RIP];
            #elif defined(__aarch64__)
                // Linux arm64 does it this way
                ip = (void*)context->uc_mcontext.pc;
            #endif
        #endif

        if (ip) {
            // We are on a platform where we can get the instruction pointer.
            *out << "Caught signal " << signalNumber << " raised at address " << ip << "; tracing with backward-cpp" << endl;
            static backward::StackTrace stack_trace;
            // With current backward-cpp we can pass the signal information and have it use the right stack.
            stack_trace.load_from(ip, 32, (void*)context, signalInfo->si_addr);
            static backward::Printer p;
            p.color_mode = backward::ColorMode::automatic;
            p.address = true;
            p.object = true;
            p.print(stack_trace, *out);

            *out << std::endl;
            *out << "Library locations:" << std::endl;

            // Now report all the objects
            for (int i = stack_trace.size(); i > 0; i--) {
                Dl_info address_library;
                if (dladdr(stack_trace[i].addr, &address_library)) {
                    *out << "#" << i << "\t";
                    report_library(*out, address_library, stack_trace[i].addr);
                    *out << std::endl;
                }
            }
        } else {
            *out << "Caught signal " << signalNumber << " at unknown address" << endl;
        }
    #endif
    
    tempStream.close();

    // Use OSC-8 to link the user to their destination.
    cerr << "ERROR: Signal "<< signalNumber << " occurred. VG has crashed. ";
    start_link(ISSUE_URL);
    cerr << "Visit ";
    cerr << ISSUE_URL;
    cerr << " to report a bug.";
    stop_link();
    cerr << endl;
    draw_br();
    dump_crash_contexts(std::cerr);
    draw_br();
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

void with_exception_handling(const std::function<void(void)>& body) {
    try {
        body();
    } catch(const std::exception& ex) {
        report_exception(ex); 
    }
}

void report_exception(const std::exception& ex) {
    #pragma omp critical (cerr)
    {
        std::cerr << std::endl;
        draw_br();
        std::cerr << "Unhandled exception of type " << typeid(ex).name() << ": " << ex.what() << std::endl;
        if (!stored_crash_context.empty()) {
            std::cerr << "Exception context: " << stored_crash_context << std::endl;
        }
    }
    abort();
}

void crash_unless_failed(const char* condition_string, const char* file, int line, const char* function) {
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
