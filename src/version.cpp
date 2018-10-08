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

// Keep the list of codenames.
// Add new codenames here
const unordered_map<string, string> Version::codenames = {
    {"v1.8.0", "Vallata"},
    {"v1.9.0", "Miglionico"},
    {"v1.10.0", "Rionero"},
    {"v1.11.0", "Cairano"}
    // Add more codenames here
};

string Version::get_version() {
    return VERSION;
}

string Version::get_release() {
    auto dash = VERSION.find('-');
    if (dash == -1) {
        // Pure tag versions have no dash
        return VERSION;
    } else {
        // Otherwise it is tag-count-ghash and the tag describes the release
        return VERSION.substr(0, dash);
    }
}

string Version::get_codename() {
    auto release = get_release();

    auto found = codenames.find(release);

    if (found == codenames.end()) {
        // No known codename for this release.
        // Return an empty string so we can just not show it.
        return "";
    } else {
        // We have a known codename!
        return found->second;
    }
}

string Version::get_short() {
    stringstream s;
    s << VERSION;

    auto codename = get_codename();
    if (!codename.empty()) {
        // Add the codename if we have one
        s << " \"" << codename << "\"";
    }

    return s.str();
}

string Version::get_long() {
    stringstream s;
    s << "vg version " << get_short() << endl;
    s << "Compiled with " << COMPILER << " on " << OS << endl;
    s << "Built by " << BUILD_USER << "@" << BUILD_HOST;
    return s.str();
}

}
