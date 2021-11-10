/**
 * \file mem_accelerator.hpp
 *
 * Implements an index for accelerating GCSA2 queries
 */

#include "mem_accelerator.hpp"
#include <sdsl/util.hpp>
#include <cmath>

namespace vg {

MEMAccelerator::MEMAccelerator(const gcsa::GCSA& gcsa_index, size_t k) : k(k)
{
    // compute the minimum width required to express the integers
    uint8_t width = 0;
    for (size_t log_cntr = 1; log_cntr <= gcsa_index.size(); log_cntr *= 2) {
        ++width;
    }
    range_table.width(max<uint8_t>(width, 1));
    // range table is initialized to size 2^(2k + 1) = 2 * 4^k
    range_table.resize(1 << (2 * k + 1));
    
    const char alphabet[5] = "ACGT";
    
    // records of (next char to query, k-mer integer encoding, range)
    vector<tuple<int64_t, int64_t, gcsa::range_type>> stack;
    stack.emplace_back(0, 0, gcsa::range_type(0, gcsa_index.size() - 1));
    
    // TODO: multithread this? probably would do it single threaded
    // to init 128 stacks or something similar
    while (!stack.empty()) {
        if (stack.size() == k + 1) {
            // we've walked the full k-mers
            range_table[2 * get<1>(stack.back())] = get<2>(stack.back()).first;
            range_table[2 * get<1>(stack.back()) + 1] = get<2>(stack.back()).second;
            stack.pop_back();
        }
        else if (get<0>(stack.back()) == 4) {
            // we've walked all the k-mers that start with this prefix
            stack.pop_back();
        }
        else {
            // extend the current range by the next character
            auto next = get<0>(stack.back())++;
            auto enc = (next << (2 * (stack.size() - 1))) | get<1>(stack.back());
            
            gcsa::range_type range;
            if (!gcsa::Range::empty(get<2>(stack.back()))) {
                range = gcsa_index.LF(get<2>(stack.back()),
                                      gcsa_index.alpha.char2comp[alphabet[next]]);
            }
            else {
                range = get<2>(stack.back());
            }
            stack.emplace_back(0, enc, range);
        }
    }
}

gcsa::range_type MEMAccelerator::memoized_LF(string::const_iterator last) const {
    int64_t enc = 0;
    for (size_t i = 0; i < k; ++i) {
        enc |= (encode(*last) << (i << 1));
        --last;
    }
    return gcsa::range_type(range_table[enc << 1], range_table[(enc << 1) | 1]);
}

}
