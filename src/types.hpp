#ifndef VG_TYPES_HPP_INCLUDED
#define VG_TYPES_HPP_INCLUDED

#include <cstdint>
#include <iostream>
#include <tuple>

#include <handlegraph/types.hpp>

/** \file
 * Contains typedefs for basic types useful for talking about graphs and
 * basic operations using them.
 */

namespace vg {

/// Represents a Node ID.
/// ID type is a 64-bit signed int.
typedef handlegraph::nid_t id_t;

/// Represents an offset along the sequence of a Node.
/// Offsets are size_t.
typedef size_t offset_t;

/// Represents an oriented position on a Node.
/// Position type: id, direction, offset.
/// Offset is counted as for as prorobuf Position, from the node's first base
/// on the forward strand, and from its last base on the reverse strand.
typedef std::tuple<id_t, bool, offset_t> pos_t;

/// Create a pos_t from a Node ID, an orientation flag, and an offset along that strand of the node.
inline pos_t make_pos_t(id_t id, bool is_rev, offset_t off) {
    return std::make_tuple(id, is_rev, off);
}

/// Extract the id of the node a pos_t is on.
inline id_t id(const pos_t& pos) {
    return std::get<0>(pos);
}

/// Return true if a pos_t is on the reverse strand of its node.
inline bool is_rev(const pos_t& pos) {
    return std::get<1>(pos);
}

/// Get the offset along the selected strand of the node from a pos_t.
inline offset_t offset(const pos_t& pos) {
    return std::get<2>(pos);
}

/// Get a reference to the Node ID of a pos_t.
inline id_t& get_id(pos_t& pos) {
    return std::get<0>(pos);
}

/// Get a reference to the reverse flag of a pos_t.
inline bool& get_is_rev(pos_t& pos) {
    return std::get<1>(pos);
}

/// Get a reference to the offset field of a pos_t, which counts along the selected strand of the node.
inline offset_t& get_offset(pos_t& pos) {
    return std::get<2>(pos);
}

/// Return true if a pos_t is unset.
inline bool is_empty(const pos_t& pos) {
    return (id(pos) == 0);
}

/// Reverse a pos_t and get a pos_t at the same **point between bases**, going the other direction.
/// To get a pos_t to the same *base*, subtract 1 from the resulting offset or call reverse_base_pos().
inline pos_t reverse(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = node_length - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

/// Reverse a pos_t and get a pos_t at the same **base**, going the other direction.
inline pos_t reverse_base_pos(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = (node_length - 1) - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

/// Print a pos_t to a stream.
inline std::ostream& operator<<(std::ostream& out, const pos_t& pos) {
    return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

} // namespace vg

#endif
