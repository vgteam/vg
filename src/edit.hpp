#ifndef VG_EDIT_HPP
#define VG_EDIT_HPP

#include "vg.pb.h"
#include <utility>
#include <iostream> // for debugging
#include "json2pb.h"

namespace vg {

using namespace std;

bool edit_is_match(const Edit& e);
bool edit_is_sub(const Edit& e);
bool edit_is_insertion(const Edit& e);
bool edit_is_deletion(const Edit& e);
pair<Edit, Edit> cut_edit_at_to(const Edit& e, size_t to_off);
pair<Edit, Edit> cut_edit_at_from(const Edit& e, size_t from_off);
// Reverse an edit and reverse complement any embedded sequence
Edit reverse_edit(const Edit& e);

}

#endif
