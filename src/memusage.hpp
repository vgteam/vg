#ifndef VG_MEMUSAGE_HPP_INCLUDED
#define VG_MEMUSAGE_HPP_INCLUDED

#include <string>

/**
 * \file memusage.hpp
 * Defines an interface to /proc/self/status and other status interfaces for debugging.
 */

namespace vg {

using namespace std;

/// Get the string value for a field in /proc/self/status by name, or "" if unsupported or not found.
string get_proc_status_value(const string& name);

/// Get the max RSS usage ever, in kb, or 0 if unsupported.
size_t get_max_rss_kb();

/// Get the max virtual memory size ever, in kb, or 0 if unsupported.
size_t get_max_vmem_kb();

/// Get the current virtual memory size, in kb, or 0 if unsupported.
size_t get_current_vmem_kb();


}

#endif
