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

bool edit_is_empty(const Edit& e) {
    return e.to_length() == 0 && e.from_length() == 0 && e.sequence().empty();
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

// breaking dependence on vg's utility.hpp
static const char complement[256] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 16
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 24
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 32
                                     'N', 'N', 'N', '$', '#', 'N', 'N', 'N', // 40 GCSA stop/start characters
                                     'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', // 48
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 56
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 64
                                     'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', // 72
                                     'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', // 80
                                     'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S', // 88
                                     'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N', // 96
                                     'N', 't', 'v', 'g', 'h', 'N', 'N', 'c', // 104
                                     'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N', // 112
                                     'N', 'q', 'y', 'w', 'a', 'a', 'b', 's', // 120
                                     'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N', // 128
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 136
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 144
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 152
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 160
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 168
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 176
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 184
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 192
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 200
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 208
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 216
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 224
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 232
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 240
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 248
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};// 256

static inline char reverse_complement(const char& c) {
    return complement[c];
}

static string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        c = complement[c];
    }
    return rc;
}


Edit reverse_complement_edit(const Edit& e) {
    // Make a reversed copy
    Edit reversed = e;
    
    // All we have to do is flip the sequence
    reversed.set_sequence(reverse_complement(e.sequence()));
    
    return reversed;
}

bool operator==(const Edit& e1, const Edit& e2) {
    return (e1.to_length() == e2.to_length())
        && (e1.from_length() == e2.from_length())
        && (e1.sequence() == e2.sequence());
}

}
