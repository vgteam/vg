#ifndef VG_LONGEST_OVERLAP_HPP_INCLUDED
#define VG_LONGEST_OVERLAP_HPP_INCLUDED

// longest_overlap.hpp: longest overlap between two strings with Z algorithm

#include <vector>
#include <algorithm>
#include <cstdint>

namespace vg {

using namespace std;

/*
 * Return the length of the longest match of a suffix of str1 with a prefix of str2
 */
template<class StringLike>
size_t longest_overlap(const StringLike& str1, const StringLike& str2);

/*
 * Return the Z array
 */
template<class StringLike>
vector<int64_t> z_algorithm(const StringLike& str);



/************
 * Template implementations
 ************/

template<class StringLike>
vector<int64_t> z_algorithm(const StringLike& str) {

    vector<int64_t> Z(str.size(), 0);

    int64_t right = 0, left = 0;
    for (int64_t i = 1; i < Z.size(); ++i) {
        if (i <= right) {
            Z[i] = min(right - i + 1, Z[i - left]);
        }
        while (i + Z[i] < Z.size() && str[Z[i]] == str[i + Z[i]]) {
            ++Z[i];
        }
        if (i + Z[i] - 1 > right) {
            left = i;
            right = i + Z[i] - 1;
        }
    }

    return Z;
}

template<class StringLike>
size_t longest_overlap(const StringLike& str1, const StringLike& str2) {

    if (str1.size() == 0 || str2.size() == 0) {
        return 0;
    }

    // string-like wrapper for concatenated string with a 0-valued sentinel separating
    class ConcatString {
    public:
        ConcatString(const StringLike& str1, const StringLike& str2) : str1(&str1), str2(&str2) {}
        size_t size() const {
            return str1->size() + str2->size() + 1;
        }
        const typename StringLike::value_type operator[](size_t i) const {
            if (i < str1->size()) {
                return (*str1)[i];
            }
            else if (i == str1->size()) {
                return typename StringLike::value_type(0);
            }
            else {
                return (*str2)[i - str1->size() - 1];
            }
        }

    private:
        const StringLike* str1 = nullptr;
        const StringLike* str2 = nullptr;
    };

    ConcatString concat(str2, str1);
    auto Z = z_algorithm(concat);

    for (size_t i = str1.size(); i != 0; --i) {
        if (Z[Z.size() - i] == i) {
            return i;
        }
    }
    return 0;
}


}

#endif
