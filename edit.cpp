#include "edit.hpp"

namespace vg {

bool edit_is_match(const Edit& e) {
    return e.from_length() == e.to_length() && e.sequence().empty();
}

bool edit_is_sub(const Edit& e) {
    return e.from_length() == e.to_length() && !e.sequence().empty();
}

bool edit_is_insertion(const Edit& e) {
    return e.from_length() == 0 && e.to_length() > 0 && !e.sequence().empty();
}

bool edit_is_deletion(const Edit& e) {
    return e.from_length() > 0 && e.to_length() == 0;
}

}
