#ifndef VG_KMP_HPP_INCLUDED
#define VG_KMP_HPP_INCLUDED

// kmp.hpp: Knuth-Morris-Pratt algorithm

#include <vector>

namespace vg {

using namespace std;

// preprocess search pattern
vector<size_t> make_prefix_suffix_table(const char* pattern, size_t len);

// return index of first match or string::npos if there is no match
size_t kmp_search(const char* text, size_t text_len,
                  const char* pattern, size_t pattern_len,
                  const vector<size_t>& prefix_suffix_table);


}

#endif
