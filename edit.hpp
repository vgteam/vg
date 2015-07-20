#ifndef VG_EDIT_HPP
#define VG_EDIT_HPP

#include "vg.pb.h"

namespace vg {

bool edit_is_match(const Edit& e);
bool edit_is_sub(const Edit& e);
bool edit_is_insertion(const Edit& e);
bool edit_is_softclip(const Edit& e);
bool edit_is_deletion(const Edit& e);
bool mapping_is_total_deletion(const Mapping& m);

}

#endif
