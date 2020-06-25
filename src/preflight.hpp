#ifndef VG_PREFLIGHT_HPP_INCLUDED
#define VG_PREFLIGHT_HPP_INCLUDED

/** \file 
 * Pre-flight checks to see if the vg binary is going to work on this system.
 * Mostly exists to check for SSE4.2 support which is still not universal.
 */

// Get standard library to identify itself by including a no-op header
#include <ciso646>

namespace vg {

#define VG_PREFLIGHT_EVERYWHERE

/// Run a preflight check to make sure that the system is usable for this build of vg.
/// Aborts with a helpful message if this is not the case.
/// We make sure to build it for a lowest-common-denominator architecture.
void preflight_check() VG_PREFLIGHT_EVERYWHERE;

}

#endif
