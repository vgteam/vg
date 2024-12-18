#ifndef VG_VERSION_HPP_INCLUDED
#define VG_VERSION_HPP_INCLUDED

// version.hpp: Version reflection information for vg builds.

#include <string>
#include <unordered_map>

namespace vg {

using namespace std;

/// Class for holding vg version and build environment information.
class Version {
public:
    /// The Git version description of this build of vg
    const static string VERSION;
    /// The compiler used to build vg
    const static string COMPILER;
    // The standard library that was used to link vg
    const static string STANDARD_LIBRARY;
    // The version of HTSlib that we saw at compile time.
    const static string HTSLIB_HEADERS;
    // The version of HTSlib that we actually linked.
    const static string HTSLIB_LIBRARY;
    /// The OS that vg was built on
    const static string OS;
    /// The user who built vg
    const static string BUILD_USER;
    /// The host that built vg
    const static string BUILD_HOST;
    
    /// Get only the version (like v1.7.0-68-g224e7625).
    static string get_version();

    /// Get the release Git tag version of vg that the current version
    /// is based on (e.g. v1.7.0-68-g224e7625 will report v1.7.0).
    static string get_release();
    
    /// Get the codename of our released version
    static string get_codename();

    /// Get a short one-line description of the current version of vg with no terminating newline.
    static string get_short();
    
    /// Get a long, potentially multi-line description of the current version
    /// of vg with no terminating newline.
    static string get_long();
private:
    // Not constructable
    Version() = delete;

    /// Store all the codenames by major version
    const static unordered_map<string, string> codenames;
};

}

#endif
