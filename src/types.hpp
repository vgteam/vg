#ifndef VG_TYPES_H
#define VG_TYPES_H

#include <tuple>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "vg.pb.h"

namespace vg {

// id type is a 64-bit signed int
typedef int64_t id_t;

// offsets are size_t
typedef size_t off_t;

// position type: id, direction, offset
typedef std::tuple<id_t, bool, off_t> pos_t;

// node to path mapping
typedef std::map<std::string, std::set<Mapping*>> NodeMapping;

//node to edges mapping
typedef std::map<id_t, std::vector<Edge*>> EdgeMapping;

}

#endif
