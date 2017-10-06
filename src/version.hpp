#ifndef VG_VERSION_HPP_INCLUDED
#define VG_VERSION_HPP_INCLUDED

// version.hpp: Store the VG Git version in its own compilation unit for fast
// rebuilds when Git state changes.

namespace vg {

using namespace std;

extern const char* VG_VERSION_STRING;

}

#endif
