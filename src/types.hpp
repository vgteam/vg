#ifndef VG_TYPES_HPP_INCLUDED
#define VG_TYPES_HPP_INCLUDED

#include <tuple>

/** \file
 * Contains typedefs for basic types useful for talking about graphs.
 */

namespace vg {

/// Represents a Node ID.
/// ID type is a 64-bit signed int.
typedef int64_t id_t;

/// Represents an offset along the sequence of a Node.
/// Offsets are size_t.
typedef size_t off_t;

/// Represents an oriented position on a Node.
/// Position type: id, direction, offset.
typedef std::tuple<id_t, bool, off_t> pos_t;

}

#endif
