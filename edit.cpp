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

pair<Edit, Edit> cut_edit_at_to(const Edit& e, size_t to_off) {
    Edit left, right;
    if (to_off > e.to_length()) {
        return make_pair(e, right);
    }
    // to-length of left portion
    size_t l = e.to_length() - to_off;
    // to-length of right portion
    size_t r = e.to_length() - l;
    if (l > e.to_length()) {
        left = e;
    } else if (edit_is_match(e)) {
        left.set_from_length(l);
        left.set_to_length(l);
        right.set_from_length(r);
        right.set_to_length(r);
    } else if (edit_is_sub(e)) {
        left.set_from_length(l);
        left.set_to_length(l);
        left.set_sequence(e.sequence().substr(0, l));
        right.set_from_length(r);
        right.set_to_length(r);
        right.set_sequence(e.sequence().substr(l));
    } else if (edit_is_insertion(e)) {
        left.set_to_length(l);
        left.set_sequence(e.sequence().substr(0, l));
        right.set_to_length(r);
        right.set_sequence(e.sequence().substr(l));
    } else if (edit_is_deletion(e)) {
        left = e;
    }
    return make_pair(left, right);
}

pair<Edit, Edit> cut_edit_at_from(const Edit& e, size_t from_off) {
    cerr << "Cutting edit " << pb2json(e) << " at " << from_off << endl;
    Edit left, right;
    if (from_off > e.from_length()) {
        return make_pair(e, right);
    }
    // from-length of left portion
    size_t l = e.from_length() - from_off;
    // from-length of right portion
    size_t r = e.from_length() - l;
    if (edit_is_match(e)) {
        left.set_from_length(l);
        left.set_to_length(l);
        right.set_from_length(r);
        right.set_to_length(r);
    } else if (edit_is_sub(e)) {
        left.set_from_length(l);
        left.set_to_length(l);
        left.set_sequence(e.sequence().substr(0, l));
        right.set_from_length(r);
        right.set_to_length(r);
        right.set_sequence(e.sequence().substr(l));
    } else if (edit_is_insertion(e)) {
        left = e;
    } else if (edit_is_deletion(e)) {
        left.set_from_length(l);
        right.set_from_length(r);
    }
    return make_pair(left, right);
}

}
