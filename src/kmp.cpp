#include "kmp.hpp"

#include <string>
#include <iostream>

namespace vg {

// TODO: these functions are very similar, could possibly merge into one?
 
vector<size_t> make_prefix_suffix_table(const char* pattern, size_t len) {
    
    vector<size_t> table(len, 0);
    
    for (size_t i = 1, j = 0; i < len;) {
        if (pattern[i] == pattern[j]) {
            ++j;
            table[i] = j;
            ++i;
        }
        else {
            if (j != 0) {
                j = table[j - 1];
            }
            else {
                table[i] = 0;
                ++i;
            }
        }
    }
    
    return table;
}

size_t kmp_search(const char* text, size_t text_len,
                  const char* pattern, size_t pattern_len,
                  const vector<size_t>& prefix_suffix_table) {
    if (text_len >= pattern_len) {
        for (size_t i = 0, j = 0, last = text_len - pattern_len; i - j <= last;) {
            if (text[i] == pattern[j]) {
                ++i;
                ++j;
                if (j == pattern_len) {
                    return i - pattern_len;
                }
            }
            else {
                if (j != 0) {
                    j = prefix_suffix_table[j - 1];
                }
                else {
                    ++i;
                }
            }
        }
    }
    return string::npos;
}


}
