#include "version.hpp"

// Get the git version macro from the build system
#include "vg_git_version.hpp"

// Do the same for the build environment info
#include "vg_environment_version.hpp"

#include <htslib/hts.h>

#include <iostream>
#include <sstream>

// If the VG_GIT_VERSION deosn't exist at all, define a placeholder
// This lets us be somewhat robust to undeterminable versions
#ifndef VG_GIT_VERSION
    #define VG_GIT_VERSION "not-from-git"
#endif

// Define a way to quote macro values.
// See https://stackoverflow.com/a/196093
#define QUOTE(arg) #arg
// We need another level to get the macro's value and not its name.
#define STR(macro) QUOTE(macro)

// Work out the standard library version
#ifdef _LIBCPP_VERSION
    #define VG_STANDARD_LIBRARY_VERSION ("libc++ " STR(_LIBCPP_VERSION))
#endif
#ifdef __GLIBCXX__
    #define VG_STANDARD_LIBRARY_VERSION ("libstd++ " STR(__GLIBCXX__))
#endif
#ifndef VG_STANDARD_LIBRARY_VERSION
    #define VG_STANDARD_LIBRARY_VERSION "unknown standard library"
#endif

// And the version of htslib
#ifndef VG_HTSLIB_VERSION
    #define VG_HTSLIB_VERSION STR(HTS_VERSION)
#endif

namespace vg {

using namespace std;

// Define all the strings as the macros' values
const string Version::VERSION = VG_GIT_VERSION;
const string Version::COMPILER = VG_COMPILER_VERSION;
const string Version::STANDARD_LIBRARY = VG_STANDARD_LIBRARY_VERSION;
const string Version::HTSLIB_HEADERS = VG_HTSLIB_VERSION;
const string Version::HTSLIB_LIBRARY(hts_version());
const string Version::OS = VG_OS;
const string Version::BUILD_USER = VG_BUILD_USER;
const string Version::BUILD_HOST = VG_BUILD_HOST;

// Keep the list of codenames.
// Add new codenames here
const unordered_map<string, string> Version::codenames = {
    {"v1.8.0", "Vallata"},
    {"v1.9.0", "Miglionico"},
    {"v1.10.0", "Rionero"},
    {"v1.11.0", "Cairano"},
    {"v1.12.0", "Parolise"},
    {"v1.12.1", "Parolise"},
    {"v1.13.0", "Moschiano"},
    {"v1.14.0", "Quadrelle"},
    {"v1.15.0", "Tufo"},
    {"v1.16.0", "Rotondi"},
    {"v1.17.0", "Candida"},
    {"v1.18.0", "Zungoli"},
    {"v1.19.0", "Tramutola"},
    {"v1.20.0", "Ginestra"},
    {"v1.21.0", "Fanano"},
    {"v1.22.0", "Rotella"},
    {"v1.23.0", "Lavello"},
    {"v1.24.0", "Montieri"},
    {"v1.25.0", "Apice"},
    {"v1.26.0", "Stornara"},
    {"v1.26.1", "Stornara"},
    {"v1.27.0", "Deliceto"},
    {"v1.27.1", "Deliceto"},
    {"v1.28.0", "Acquafredda"},
    {"v1.29.0", "Sospiro"},
    {"v1.30.0", "Carentino"},
    {"v1.31.0", "Caffaraccia"},
    {"v1.32.0", "Sedlo"},
    {"v1.33.0", "Moscona"},
    {"v1.34.0", "Arguello"},
    {"v1.35.0", "Ghizzano"},
    {"v1.36.0", "Cibottola"},
    {"v1.37.0", "Monchio"},
    {"v1.38.0", "Canossa"},
    {"v1.39.0", "Runzi"},
    {"v1.40.0", "Suardi"},
    {"v1.41.0", "Salmour"},
    {"v1.42.0", "Obolo"},
    {"v1.43.0", "Barisano"},
    {"v1.44.0", "Solara"},
    {"v1.45.0", "Alpicella"},
    {"v1.46.0", "Altamura"},
    {"v1.47.0", "Ostuni"},
    {"v1.48.0", "Gallipoli"},
    {"v1.49.0", "Peschici"},
    {"v1.50.0", "Monopoli"},
    {"v1.50.1", "Monopoli"},
    {"v1.51.0", "Quellenhof"},
    {"v1.52.0", "Bozen"},
    {"v1.53.0", "Valmontone"},
    {"v1.54.0", "Parafada"},
    {"v1.55.0", "Bernolda"},
    {"v1.56.0", "Collalto"},
    {"v1.57.0", "Franchini"},
    {"v1.58.0", "Cartari"},
    {"v1.59.0", "Casatico"},
    {"v1.60.0", "Annicco"},
    {"v1.61.0", "Plodio"},
    {"v1.62.0", "Ranzano"},
    {"v1.63.0", "Boccaleone"},
    {"v1.63.1", "Boccaleone"},
    {"v1.64.0", "Vibbiana"},
    {"v1.64.1", "Vibbiana"},
    {"v1.65.0", "Carfon"},
    {"v1.65.1", "Carfon"},
    {"v1.66.0", "Navetta"},
    {"v1.67.0", "Vetria"},
    {"v1.68.0", "Rimbocchi"},
    {"v1.69.0", "Bologna"},
    {"v1.70.0", "Zebedassi"}
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
    s << "Linked against " << STANDARD_LIBRARY << endl;
    s << "Using HTSlib headers " << HTSLIB_HEADERS << ", library " << HTSLIB_LIBRARY << endl;
    s << "Built by " << BUILD_USER << "@" << BUILD_HOST;
    return s.str();
}

}
