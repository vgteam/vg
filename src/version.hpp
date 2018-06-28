#ifndef VG_VERSION_HPP_INCLUDED
#define VG_VERSION_HPP_INCLUDED

// version.hpp: Version reflection information for vg builds.

#include <string>

namespace vg {

using namespace std;

/// Class for holding vg version and build environment information.
class Version {
public:
    /// The Git version description of this build of vg
    const static string VERSION;
    /// The compiler used to build vg
    const static string COMPILER;
    /// The OS that vg was built on
    const static string OS;
    /// The user who built vg
    const static string BUILD_USER;
    /// The host that built vg
    const static string BUILD_HOST;
    
    /// Get a short description of the current version of vg.
    static const string& get_short();
    
    /// Get a long, potentially multi-line description of the current version
    /// of vg with no terminating newline.
    static string get_long();
private:
    // Not constructable
    Version() = delete;
};

}

#endif
