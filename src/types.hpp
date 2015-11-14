#ifndef VG_TYPES_H
#define VG_TYPES_H

#include <tuple>

namespace vg {

// id type is a 64-bit signed int
typedef int64_t id_t;

// offsets are size_t
typedef size_t off_t;

// position type: id, direction, offset
typedef std::tuple<id_t, bool, off_t> pos_t;

}

#endif
