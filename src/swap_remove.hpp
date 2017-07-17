#ifndef VG_SWAP_REMOVE_HPP_INCLUDED
#define VG_SWAP_REMOVE_HPP_INCLUDED

#include <vector>
#include <algorithm>

template<typename T>
bool swap_remove(std::vector<T>& v, const T& e) {
    bool swapped = false;
    for (typename std::vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
        if (e == *i) {
            std::swap(*i, v.back());
            swapped = true;
            break;
        }
    }
    if (swapped) {
        v.pop_back();
        return true;
    } else {
        return false;
    }
}

#endif
