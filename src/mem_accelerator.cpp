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
    // Compute the minimum width required to express the integers.
    // We need to be able to represent numbers up to a bit over the index size,
    // for some reason. We're not sure how much over the index size is needed.
    // Start with a more or less arbitrary guess of the max needed value
    size_t max_needed_value_guess = gcsa_index.size() + 100;
    
    // Compute the bit width needed to represent it
    uint8_t width = 0;
    // And the number too big to be represented
    size_t too_big = 1;
    while (too_big <= max_needed_value_guess) {
        ++width;
        too_big *= 2;
    }
    
    // We use this accessor to set a value and expand the vector width if we need to.
    auto set_range_table = [&](size_t offset, int64_t value) {
        if (too_big <= value) {
            while (too_big <= value) {
                ++width;
                too_big *= 2;
            }
            // This is weird; we should sort this out when we work out exactly
            // what the limits on these range values really are.
            std::cerr << "warning [vg::MEMAccelerator]: expanding vector width to hold value " << value << " for GCSA index size " << gcsa_index.size() << std::endl;
            sdsl::util::expand_width(range_table, width);
        }
        range_table[offset] = value;
    };
    
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
            set_range_table(2 * get<1>(stack.back()), get<2>(stack.back()).first);
            set_range_table(2 * get<1>(stack.back()) + 1, get<2>(stack.back()).second);
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
