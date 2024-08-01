#ifndef VG_MINWIDTH_INT_HPP_INCLUDED
#define VG_MINWIDTH_INT_HPP_INCLUDED

#include <vector>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <limits>



/** \file min_width_int_vector.hpp
 * Methods for storing a vector of integers with minimal bit width
 */

namespace vg{
using namespace std;

/* A struct to store a vector of integers with minimal bit width
 */
struct min_width_int_vector_t {

    private:

    /// How many bits are used to store the bit width used
    /// This is needed for serializing
    const static size_t BIT_WIDTH_WIDTH = 8;

    /// The bit width that is being used to store the integers
    uint8_t width;

    ///The actual data stored in the vector
    std::vector<bool> data;

    public:

    min_width_int_vector_t () {
        width = 0;
    }

    min_width_int_vector_t (size_t w) {
        width = w;
    }


    ///Make this a copy of input_data 
    ///If maxval is set, then this is the maximum value in the input data, 
    /// or the maximum value to be stored with the bitwidth
    ///If there is no max_val and the width has not already been set, get the
    /// width from the maximum value in input_data
    void from_vector(const vector<size_t>& input_data, size_t max_val = 0);

    ///Add a value to the end of the vector
    void push_back(size_t val);

    ///How long is the vector
    size_t size() const;

    ///Get the value at the given index
    size_t at(size_t index) const;

    ///Check what the bit width is
    // This is a size_t because it's blank when I try to write it to stderr
    size_t get_bit_width() const { return (size_t) width;}

    ///How many bits are we using total
    size_t get_bit_count() const { return data.size(); }

    ///////////Access the bit vector itself for serializing
    bool bit_at(size_t i) const {return data[i];}
    void set_bitvector_length(size_t l) {data.resize(l);}
    void set_bit_at(size_t i) {data[i] = true;}
    void set_bit_width(size_t w) {width = w;}


    ///Equality operator
    //TODO: This isn't actually checking the values- the widths could be different but still represent the same vectors.
    //      but that would be pretty slow to check so leave it
    inline bool operator==(const min_width_int_vector_t& other) const {
        return width == other.width && data == other.data;
    }

};
}
#endif
