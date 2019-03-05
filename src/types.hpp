#ifndef VG_TYPES_HPP_INCLUDED
#define VG_TYPES_HPP_INCLUDED

#include <cstdint>
#include <tuple>

#include <handlegraph/types.hpp>

/** \file
 * Contains typedefs for basic types useful for talking about graphs.
 */

namespace vg {

/// Represents a Node ID.
/// ID type is a 64-bit signed int.
typedef handlegraph::nid_t id_t;

/// Represents an offset along the sequence of a Node.
/// Offsets are size_t.
typedef size_t off_t;

/// Represents an oriented position on a Node.
/// Position type: id, direction, offset.
/// Offset is counted as for as prorobuf Position, from the node's first base
/// on the forward strand, and from its last base on the reverse strand.
typedef std::tuple<id_t, bool, off_t> pos_t;

}

#endif
