#include "preflight.hpp"

#include <cstdio>
#include <cstdlib>

#ifdef __x86_64__
#include <cpuid.h>
#endif

namespace vg {

using namespace std;

// Define our error message so we can be clear that we're working with its address.
const static char* PREFLIGHT_FAIL_MESSAGE = "error[vg::preflight_check]: The CPU does not support SSE4.2 instructions. VG cannot run here. Please use a system with SSE4.2 support.\n";
// Define function pointers for the standard library functions we are going to
// call to be super sure they aren't inlined. If they are marked always-inline
// but also not marked to be built for the arch of the calling function we get
// a "target specific option mismatch" error.
static auto* fputs_ptr = &fputs;
static auto* exit_ptr = &exit;

void preflight_check() {
    // This whole function needs to run without nice things like C++ allocators
    // or std::endl, which are likely to be both always-inline and not compiled
    // for the lowest common denominator architecture.
    //
    // TODO: Build the whole compilation unit for the lowest common denominator
    // architecture?

    bool arch_ok = true;

#ifdef __x86_64__
    // We assume we are on x86_64 on POSIX (and not Windows).
    // We use the method of dlib's dlib/simd/simd_check.h

    // Define a place to put the cpuid info
    unsigned int cpuid_info[4];

    // Call cpuid function 1 (which reports SSE4.2, and other stuff up to original AVX)
    __cpuid(1, cpuid_info[0], cpuid_info[1], cpuid_info[2], cpuid_info[3]);

    // Bit 20 of result 2 is the SSE 4.2 flag.
    bool have_sse4_2 = cpuid_info[2] & (1 << 20);

    arch_ok &= have_sse4_2;
#endif
    // If not on x86_64, we are probably on ARM and using fake SSE anyway.

    if (!arch_ok) {
        // Call the function addresses with normal call instructions.
        // Hope we didn't statically link libc, or that they work here.
        (*fputs_ptr)(PREFLIGHT_FAIL_MESSAGE, stderr);
        (*exit_ptr)(1);
    }
}

}
