#ifndef VG_VARINT_HPP_INCLUDED
#define VG_VARINT_HPP_INCLUDED

#include <vector>
#include <utility>
#include <cstdint>

/** \file varint.hpp
 * Methods for storing a vector of integers with variable bit width
 * Implements protobuf's varints
 */
#define DEBUG_VARINT

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
        const inline std::pair<size_t, size_t> get_value_and_next_index(size_t index);
    
        private:
        //The actual data stored in the vector
        std::vector<uint8_t> data;

        const static size_t USABLE_BITS = 7;
        //01111111
        const static uint8_t MAX_VALUE = (1 << USABLE_BITS) - 1;
    
    };

void write_byte_as_bits_to_stderr(size_t value) {
    cerr << ((value & (1<<7)) ? "1" : "0") 
         << ((value & (1<<6)) ? "1" : "0") 
         << ((value & (1<<5)) ? "1" : "0") 
         << ((value & (1<<4)) ? "1" : "0") 
         << ((value & (1<<3)) ? "1" : "0") 
         << ((value & (1<<2)) ? "1" : "0") 
         << ((value & (1<<1)) ? "1" : "0") 
         << ((value & (1<<0)) ? "1" : "0");
}

    /*The values get stored in chunks of 7 bits, with the 7 least significant bits first.
     * The first bit in each byte of the vector data indicates whether the next byte is part
     * of the same value (1 to continue, 0 if the current byte is the last in the integer)
     * TODO: This assumes that everything is big-endian, which may not be true?
     */
    
    void varint_vector_t::add_value(size_t value) {
        if (value == 0) {
            //If the value is 0, then the 0 tag to end the integer and 0 for the value 
#ifdef DEBUG_VARINT
                cerr <<"adding " << data.size() << ": 0" << endl;
#endif
            data.push_back(0);
            return;
        }
        while (value != 0) {
            if (value < MAX_VALUE) {
                //If the remainder of the integer can be stored in 7 bits
                //then it gets stored with a 0 as the first bit
#ifdef DEBUG_VARINT
                cerr <<"adding " << data.size() << ": ";
                write_byte_as_bits_to_stderr(value);
                cerr << endl;
#endif
                data.push_back(value);
            } else {
                //Otherwise, store a byte with a 1 as the first bit, and then the 
                //7 least significant bits of value
#ifdef DEBUG_VARINT
                cerr << "adding " << data.size() << ": ";
                write_byte_as_bits_to_stderr((1<<USABLE_BITS) | ( MAX_VALUE & value));
                cerr << endl;
#endif
                data.push_back((1<<USABLE_BITS) | ( MAX_VALUE & value));
            }

            //right shift value to get rid of the last 7 bits 
            value = value >> USABLE_BITS;
        }
        
        return;
    }
    
    //TODO: What to do if its empty?
    const inline std::pair<size_t, size_t> varint_vector_t::get_value_and_next_index(size_t index) {
        if (index >= data.size()) {
            throw runtime_error("Accessing value past the end of a varint vector");
        }

        //Value to return
        size_t value = 0;
        //How many chunks have we seen so far
        size_t chunk_count = 0;

        //TODO: Shouldn't have to check the size of the array because the last thing should have a 0 in front of it anyway
        while (index < (data.size()-1) && (data[index]>>USABLE_BITS) == 1) {
#ifdef DEBUG_VARINT
            cerr << "retrieving: " << index << ": ";
            write_byte_as_bits_to_stderr(data[index]);
            cerr << endl;
#endif
            //For each chunk, add the 7 bits from the current index to value
            //TODO: I'd like to not have to explicitly make a new size_t but reinterpret_cast doesn't compile and it'll cut off after 32 bits otherwise
            size_t to_add = (data[index] & MAX_VALUE);
            value |= (to_add << (USABLE_BITS*chunk_count));

            //Increment the current index and the number of things we've added
            index++;
            chunk_count++;
        }

        //After the loop, either the index points to the last thing or the current byte that index
        //points to starts with a 0, indicating that it's the last chunk of the current value
#ifdef DEBUG_VARINT
        cerr << "retrieving: " << index << ": ";
        write_byte_as_bits_to_stderr(data[index]);
        cerr << endl;
        write_byte_as_bits_to_stderr((data[index] & MAX_VALUE));
        cerr << " " << (USABLE_BITS*chunk_count) << endl;
#endif
            size_t to_add = (data[index] & MAX_VALUE);
            value |= (to_add << (USABLE_BITS*chunk_count));

        index++;

        return std::make_pair(value, index);
    }

}
#endif
