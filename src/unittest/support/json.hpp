#ifndef VG_UNITTEST_JSON_HPP_INCLUDED
#define VG_UNITTEST_JSON_HPP_INCLUDED
/** \file json.hpp
 * Utilities for working with JSON data in test cases.
 */

#include "handle.hpp"
#include <utility>


namespace vg {
namespace unittest {

/// Create a handlegraph from vg Protobuf JSON.
std::unique_ptr<MutablePathMutableHandleGraph> json_to_graph(const std::string& json);


}
}

#endif
