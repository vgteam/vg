#include "version.hpp"

// Get the git version macro from the build system
#include "vg_git_version.hpp"

// Do the same for the build environment info
#include "vg_environment_version.hpp"

#include <iostream>
#include <sstream>

namespace vg {

using namespace std;

// Define all the strings as the macros' values
const string Version::VERSION = VG_GIT_VERSION;
const string Version::COMPILER = VG_COMPILER_VERSION;
const string Version::OS = VG_OS;
const string Version::BUILD_USER = VG_BUILD_USER;
const string Version::BUILD_HOST = VG_BUILD_HOST;

const string& Version::get_short() {
    return VERSION;
}

string Version::get_long() {
    stringstream s;
    s << "vg version " << VERSION << endl;
    s << "Compiled with " << COMPILER << " on " << OS << endl;
    s << "Built by " << BUILD_USER << "@" << BUILD_HOST;
    return s.str();
}

}
