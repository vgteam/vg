#include "min_width_int_vector.hpp"
#include <iostream>
#include <sstream>
#include <limits>

//#define DEBUG_MININT

namespace vg {
using namespace std;

void min_width_int_vector_t::from_vector(const vector<size_t>& input_data, size_t max_val) {
    if (max_val != 0) {
        width = std::max(width, 1 + (size_t)std::floor(std::log2(max_val)));
    } else if (width == 0) {
        //If we haven't already set the width, find it from the max value of the input data
        for (const size_t& x : input_data) {
            max_val = std::max(x, max_val);
        }
        width = 1 + (size_t)std::floor(std::log2(max_val));
    }
    data.reserve(input_data.size()*width);

    for (const size_t& x : input_data) {
        push_back(x);
    }
}

void min_width_int_vector_t::push_back(size_t val) {
#ifdef DEBUG_MININT
    assert(width >= 1 + (size_t)std::floor(std::log2(val)));
#endif
    for (size_t i = 0 ; i < width ; i++) {
        data.emplace_back(val & (1 << (width - i - 1)));
    }                
                
}

size_t min_width_int_vector_t::size() const {
    return data.size() / width;
}
size_t min_width_int_vector_t::at(size_t index) const {
    size_t result = 0;
    size_t start_index = index * width;
    for (size_t i = 0 ; i < width ; i++) {
        if (data[i + start_index]) {
            result |= (1 << (width - i - 1));
        }
    }
    return result; 
}


}
