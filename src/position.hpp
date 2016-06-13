#ifndef VG_POS_H
#define VG_POS_H

#include "vg.pb.h"
#include "types.hpp"
#include <iostream>

namespace vg {

using namespace std;

id_t id(const pos_t& pos);
bool is_rev(const pos_t& pos);
off_t offset(const pos_t& pos);
id_t& get_id(pos_t& pos);
bool& get_is_rev(pos_t& pos);
off_t& get_offset(pos_t& pos);
pos_t reverse(const pos_t& pos, size_t node_length);
Position reverse(const Position& pos, size_t node_length);
ostream& operator<<(ostream& out, const pos_t& pos);
pos_t make_pos_t(const Position& pos);
pos_t make_pos_t(id_t id, bool is_rev, off_t off);
Position make_position(const pos_t& pos);
Position make_position(id_t id, bool is_rev, off_t off);

}

#endif
