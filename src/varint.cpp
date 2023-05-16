#include "varint.hpp"
#include <iostream>
#include <sstream>
#include <limits>

//#define DEBUG_VARINT

namespace vg {
using namespace std;

#ifdef DEBUG_VARINT
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
#endif

/*The values get stored in chunks of 7 bits, with the 7 least significant bits first.
 * The first bit in each byte of the vector data indicates whether the next byte is part
 * of the same value (1 to continue, 0 if the current byte is the last in the integer)
 * TODO: This assumes that everything is big-endian, which may not be true?
 */

void varint_vector_t::add_value(size_t value) {
#ifdef DEBUG_VARINT
    cerr << "Set varint_vector(" << (void*)this << ")[" << data.size() << "] = " << value << endl;
#endif
    if (value == 0) {
        //If the value is 0, then the 0 tag to end the integer and 0 for the value 
#ifdef DEBUG_VARINT
            cerr <<"adding " << data.size() << ": 0" << endl;
#endif
        data.push_back(0);
        return;
    }
    while (value != 0) {
        if (value <= MAX_VALUE) {
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
std::pair<size_t, size_t> varint_vector_t::get_value_and_next_index(size_t index) const {
#ifdef DEBUG_VARINT
    size_t original_index = index;
#endif
    if (index >= data.size()) {
        std::stringstream ss;
        // Note that this is the address of the varint_vector_t and not its data.
        ss << "Accessing value at " << index << " past the end of a varint vector size " << data.size() << " at " << (void*) this;
        throw runtime_error(ss.str());
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

    //If this was the last thing in the list, then return std::numeric_limits<size_t>::max() as
    //the next index
    if (index == data.size()) {
        index = std::numeric_limits<size_t>::max();
    }

#ifdef DEBUG_VARINT
    cerr << "Found varint_vector(" << (void*)this << ")[" << original_index << "] = " << value << ", " << index << endl;
#endif

    return std::make_pair(value, index);
}

void varint_vector_t::print_self() const {
    for (const auto& byte : data) {
        cerr << (static_cast<size_t>(byte)) <<  ":  "
             << ((byte & (1<<7)) ? "1" : "0") 
             << ((byte & (1<<6)) ? "1" : "0") 
             << ((byte & (1<<5)) ? "1" : "0") 
             << ((byte & (1<<4)) ? "1" : "0") 
             << ((byte & (1<<3)) ? "1" : "0") 
             << ((byte & (1<<2)) ? "1" : "0") 
             << ((byte & (1<<1)) ? "1" : "0") 
             << ((byte & (1<<0)) ? "1" : "0") << endl;
    }
}

std::vector<size_t> varint_vector_t::to_vector() const {
    std::vector<size_t> to_return;

    std::pair<size_t, size_t> value_and_index = {0, 0};

    while (value_and_index.second < data.size()) {
        // Until we hit the end of our data, decode values and store them.
        value_and_index = get_value_and_next_index(value_and_index.second);
        to_return.push_back(value_and_index.first);
    }

    return to_return;
}

void varint_vector_t::from_vector(const std::vector<size_t>& values) {
    // Throw away anything we have already
    data.clear();
    for (auto& v : values) {
        // And encode all the values we were given
        add_value(v);
    }
}

}
