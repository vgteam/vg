#ifndef VG_VARINT_HPP_INCLUDED
#define VG_VARINT_HPP_INCLUDED

#include <vector>
#include <utility>
#include <cstdint>

/** \file varint.hpp
 * Methods for storing a vector of integers with variable bit width
 * Implements protobuf's varints
 */

namespace vg{
using namespace std;

    /* A struct to store a vector of integers with variable bit width
     * Values can only be accessed in order, and only added to the end of the vector
     */
    struct varint_vector_t {
        public:
    
        //Add an integer value to the end of the varint vector
        void add_value(size_t value);
    
        //Get the integer at the given index. 
        //Index refers to the index in the vector of bytes, not the nth value stored in the vector
        //Also return the index of the next value
        //Returns std::numeric_limits<size_t>::max() as the next index if the current index is the 
        //last thing in the vector
        std::pair<size_t, size_t> get_value_and_next_index(size_t index) const;

        ///Equality operator
        inline bool operator== (const varint_vector_t& other ) const{
            return data == other.data;
        }
    
        private:
        //The actual data stored in the vector
        std::vector<uint8_t> data;

        const static size_t USABLE_BITS = 7;
        //01111111
        const static uint8_t MAX_VALUE = (1 << USABLE_BITS) - 1;


    
    };
}
#endif
