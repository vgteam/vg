#include "min_width_int_vector.hpp"

//#define DEBUG_MININT

namespace vg {
using namespace std;

void min_width_int_vector_t::from_vector(const vector<size_t>& input_data, size_t max_val) {
#ifdef DEBUG_MININT
    cerr << "get minint vector from int vector " << endl;
#endif
    if (max_val != 0) {
#ifdef DEBUG_MININT
        cerr << "Get width from max value " << max_val << " bigger of " << ((size_t) width) << " and " << (std::floor(std::log2(max_val)) + 1) << endl;
#endif
        width = (uint8_t) std::max((size_t) width, (size_t)(std::floor(std::log2((float) max_val)) + 1));
    } else if (width == 0) {
        //If we haven't already set the width, find it from the max value of the input data
        for (const size_t& x : input_data) {
            max_val = std::max(x, max_val);
        }
#ifdef DEBUG_MININT
        cerr << "Found max value " << max_val << " and got width " << width << endl;
#endif
        width = 1 + (size_t)std::floor(std::log2((float) max_val));
    }
#ifdef DEBUG_MININT
    for (size_t x : input_data) {
        cerr << x << " ";
    }
    for (size_t x : input_data) {
        assert( width >= (uint8_t)(std::floor(std::log2(x)) + 1));
    }
#endif
    data.reserve(input_data.size()*width);

    for (const size_t& x : input_data) {
        push_back(x);
    }
}



void min_width_int_vector_t::push_back(size_t val) {
#ifdef DEBUG_MININT
    assert(width >= (uint8_t) (1 + (size_t)std::floor(std::log2(val))));
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
