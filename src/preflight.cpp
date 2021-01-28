#include "preflight.hpp"

#include <iostream>
#include <cstdlib>

#ifdef __x86_64__
#include <cpuid.h>
#endif

namespace vg {

using namespace std;

void preflight_check() {
    
#ifdef __x86_64__
    // We assume we are on x86_64 on POSIX (and not Windows).
    // We use the method of dlib's dlib/simd/simd_check.h
    
    // Define a place to put the cpuid info
    unsigned int cpuid_info[4];
    
    // Call cpuid function 1 (which reports SSE4.2, and other stuff up to original AVX)
    __cpuid(1, cpuid_info[0], cpuid_info[1], cpuid_info[2], cpuid_info[3]);
    
    // Bit 20 of result 2 is the SSE 4.2 flag.
    bool have_sse4_2 = cpuid_info[2] & (1 << 20);
    
    if (!have_sse4_2) {
        cerr << "error[vg::preflight_check]: The CPU does not support SSE4.2 instructions. VG cannot run here. "
            << "Please use a system with SSE4.2 support." << endl;
        exit(1);
    }
#endif
    // If not on x86_64, we are probably on ARM and using fake SSE anyway.
    
}

}
