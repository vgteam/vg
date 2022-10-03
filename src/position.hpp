#ifndef VG_POSITION_HPP_INCLUDED
#define VG_POSITION_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "utility.hpp"
#include "vg/io/json2pb.h"
#include <gcsa/gcsa.h>
#include "handle.hpp"

/** \file 
 * Functions for working with Positions and `pos_t`s.
 */

namespace vg {

using namespace std;



// duplicative with pos_t, but this will make it so much easier
// to refactor -- can always eliminate later
class position_t {
public:
    position_t() = default;
    position_t(const position_t&) = default;
    position_t(position_t&&) = default;
    ~position_t() = default;
    position_t& operator=(const position_t&) = default;
    position_t& operator=(position_t&&) = default;
    inline int64_t node_id() const;
    inline void set_node_id(int64_t i);
    inline int64_t offset() const;
    inline void set_offset(int64_t o);
    inline bool is_reverse() const;
    inline void set_is_reverse(bool r);
    inline bool operator==(const position_t& other) const;
    inline bool operator!=(const position_t& other) const;
private:
    int64_t _node_id;
    int64_t _offset;
    bool _is_reverse;
};

/// Reverse a Position and get a Position at the same **point between bases**, going the other direction.
/// To get a Position to the same *base*, subtract 1 from the resulting offset.
Position reverse(const Position& pos, size_t node_length);
/// Convert a Position to a (much smaller) pos_t.
pos_t make_pos_t(const Position& pos);
pos_t make_pos_t(const position_t& pos);
/// Create a pos_t from a gcsa node
pos_t make_pos_t(gcsa::node_type node);
/// Convert a pos_t to a Position.
Position make_position(const pos_t& pos);
/// Create a Position from a Node ID, an orientation flag, and an offset along that strand of the node.
Position make_position(id_t id, bool is_rev, vg::off_t off);
/// Make a Position from a gcsa node
Position make_position(gcsa::node_type node);
/// Find the min distance in the path offsets where the path orientation is the same and different
pair<int64_t, int64_t> min_oriented_distances(const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& path_offsets1,
                                              const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& path_offsets2);

string debug_string(const position_t& pos);
void from_proto_position(const Position& from, position_t& to);

/*
 * position_t
 */
inline int64_t position_t::node_id() const {
    return _node_id;
}
inline void position_t::set_node_id(int64_t i) {
    _node_id = i;
}
inline int64_t position_t::offset() const {
    return _offset;
}
inline void position_t::set_offset(int64_t o) {
    _offset = o;
}
inline bool position_t::is_reverse() const {
    return _is_reverse;
}
inline void position_t::set_is_reverse(bool r) {
    _is_reverse = r;
}

inline bool position_t::operator==(const position_t& other) const {
    return (_node_id == other._node_id
            && _is_reverse == other._is_reverse
            && _offset == other._offset);
}

inline bool position_t::operator!=(const position_t& other) const {
    return !(*this == other);
}
}

#endif
