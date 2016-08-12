#include "version.hpp"

// Get the actual macro from the build system
#include "vg_git_version.hpp"

// Make sure the version macro is a thing
#ifndef VG_GIT_VERSION
    #define VG_GIT_VERSION "missing"
#endif

namespace vg {

using namespace std;

// Define the constant equal to the macro.
const char* VG_VERSION_STRING = VG_GIT_VERSION;

}
