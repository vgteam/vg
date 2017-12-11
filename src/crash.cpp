#include "crash.hpp"

// Needed for automatic name demangling, but not all that portable
#include <cxxabi.h>
#include <execinfo.h>
// Needed to hack contexts for signal traces
#ifndef _XOPEN_SOURCE
    #define _XOPEN_SOURCE
    #include <ucontext.h>
    #undef _XOPEN_SOURCE
#else
    #include <ucontext.h>
#endif
// Needed to generate stacktraces ourselves
#include <dlfcn.h>

// iostream wants this on Travis on Mac
#include <pthread.h>

// We need strcmp
#include <cstring>

#include <cstdlib>
#include <string>
#include <iostream>

namespace vg {

// Demangle the name in thsi stack trace frame if we can find the API to do so.
string demangle_frame(string mangled) {
    // Demangle the name in a stack trace line as seen at
    // <http://panthema.net/2008/0901-stacktrace-demangled/>
    
    // Name is module(function+offset) [address] in standard format.
    // For example: 
    // ../createIndex/createIndex(_Z12make_tempdirv+0x1a4) [0x46e8f4]
    // We need to find the start and end of the function part. Sometimes
    // the parens can just be empty, so we need to handle that too.
    
    // Where is the close paren, reading from the right?
    size_t closeParen = 0;
    // Where is the plus, reading from the right? If there is no plus in
    // the parens, set to 0.
    size_t plus = 0;
    // Where is the open paren, reading from right to left?
    size_t openParen = 0;
    
    for(size_t j = mangled.size() - 1; j != (size_t) -1; j--) {
        // Scan from right to left.
        
        if(closeParen == 0 && mangled[j] == ')') {
            // We found the rightmost close paren
            closeParen = j;
        } else if(j < closeParen && plus == 0 && mangled[j] == '+') {
            // We found the + to the left of the close paren.
            plus = j;
        } else if(j < closeParen && openParen == 0 && mangled[j] == '(') {
            // We found the open paren to the left of the close paren.
            openParen = j;
            
            // We're done parsing.
            break;
        }
    }
    
    if(openParen == 0 || closeParen == 0 || plus == 0) {
        // We couldn't pull out a name and address. Either we have a
        // nonstandard format or we have empty parens.
        
        // Just use the default trace message
        return mangled;
    } else {
        // We did parse out stuff!
        
        // Take everything before the open paren.
        string demangled = mangled.substr(0, openParen + 1);
        
        // Grab the function name
        string functionName = mangled.substr(openParen + 1, plus - (openParen + 1));
        
        // Make a place for the demangling function to save its status
        int status;
        
        // Do the demangling
        char* demangledName = abi::__cxa_demangle(functionName.c_str(), NULL, NULL, &status);
        
        if(status != 0) {
            // If we couldn't demangle the name, just use the mangled name.
            return mangled;
        }
        
        // Add the (probably) demangled name, a "+", and the rest of the
        // message.
        demangled += string(demangledName) + "+" + mangled.substr(plus + 1);
        
        if(status == 0) {
            // We got a demangled name we need to clean up.
            free(demangledName);
        }
        
        return demangled;
    }
}

void stacktrace_with_backtrace_and_exit(int signalNumber) {
    // How many frames can we handle?
    const size_t MAX_FRAMES = 100;
    
    // This holds the stack frames
    void *frames[MAX_FRAMES];
    
    // And this holds how many there actually are, which comes out of the
    // function that gets the frames.
    size_t framesUsed = backtrace(frames, MAX_FRAMES);
    
    cerr << "Stack trace from backtrace() for signal " << signalNumber << ":" << endl;
        
    char** traceMessages = backtrace_symbols(frames, framesUsed);
    
    for(size_t i = 0; i < framesUsed; i++) {
        // Print a demangled version of every frame            
        cerr << demangle_frame(traceMessages[i]) << endl;
        // Separate frames because damangled can be long.
        cerr << "=================" << endl;
    }
    
    // Free our stacktrace memory.
    free(traceMessages);
    
    exit(signalNumber + 128);
}

void emit_stacktrace(int signalNumber, siginfo_t *signalInfo, void *signalContext) {
    // This holds the context that the signal came from, including registers and stuff
    ucontext_t* context = (ucontext_t*) signalContext;
    
    cerr << "Got signal " << signalNumber << endl;
    
    cerr << "Manual stack trace:" << endl;
    
    #ifdef __APPLE__
        // On OS X we need to do our own stack tracing since backtrace() doesn't seem to work
        
        // Now we compute our own stack trace, because backtrace() isn't so good on OS X.
        // We operate on the same principles as <https://stackoverflow.com/a/5426269>
        
        // TODO: This assumes x86
        // Fetch out the registers
        // We model IP as a pointer to void (i.e. into code)
        void* ip;
        // We model BP as an array of two things: previous BP, and previous IP.
        void** bp;
        
        // OS X 64 bit does it this way
        ip = (void*)context->uc_mcontext->__ss.__rip;
        bp = (void**)context->uc_mcontext->__ss.__rbp;
        
        // If this were Linux it would be (void*)context->uc_mcontext.gregs[REG_RIP] or REG_RBP
        
        // Allocate a place to keep the dynamic library info for the address the stack is executing at
        Dl_info address_library;
        
        cerr << endl;
        cerr << "Next ip: " << ip << " Next bp: " << bp << endl;
        while (true) {
            // Work out where the ip is.
            if (!dladdr(ip, &address_library)) {
                // This address doesn't belong to anything!
                cerr << "Stack leaves code at ip=" << ip << endl;
                break;
            }

            if (address_library.dli_sname != nullptr) {
                // Make a place for the demangling function to save its status
                int status;
                
                // Do the demangling
                char* demangledName = abi::__cxa_demangle(address_library.dli_sname, NULL, NULL, &status);
                
                if(status == 0) {
                    // Successfully demangled
                    cerr << "Address " << ip << " in demangled symbol " << demangledName
                        << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_saddr))
                        << ", in library " << address_library.dli_fname
                        << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_fbase)) << endl;
                    free(demangledName);
                } else {
                    // Leave mangled
                    cerr << "Address " << ip << " in mangled symbol " << address_library.dli_sname
                        << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_saddr))
                        << ", in library " << address_library.dli_fname
                        << " at offset " << (void*)((size_t)ip - ((size_t)address_library.dli_fbase)) << endl;
                }
                
                #ifdef VG_DO_ATOS
                    // Try running atos to print a line number. This can be slow so we don't do it by default.
                    stringstream command;
                    
                    command << "atos -o " << address_library.dli_fname << " -l " << address_library.dli_fbase << " " << ip;
                    cerr << "Running " << command.str() << "..." << endl;
                    system(command.str().c_str());
                #endif
                
            } else {
                cerr << "Address " << ip << " out of symbol in library " << address_library.dli_fname << endl;
            }

            if(address_library.dli_sname != nullptr && !strcmp(address_library.dli_sname, "main")) {
                cerr << "Stack hit main" << endl;
                break;
            }

            if (bp != nullptr) {
                // Simulate a return
                ip = bp[1];
                bp = (void**) bp[0];
            } else {
                break;
            }
            
            cerr << endl;
            cerr << "Next ip: " << ip << " Next bp: " << bp << endl;
        }
        
        cerr << "Stack trace complete" << endl;
        cerr << endl;
        
        // Make sure to exit with the right code
        exit(signalNumber + 128);
        
    #else
        // On Linux we should just use Backtrace
        stacktrace_with_backtrace_and_exit(signalNumber);
    #endif
    
    
    
    
}

}
