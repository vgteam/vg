#ifndef VG_ENDIANNESS_HPP_INCLUDED
#define VG_ENDIANNESS_HPP_INCLUDED

/** \file endianness.hpp
 *  Methods for converting endianness in integers
 */

#include <cstdlib>
#include <type_traits>

namespace vg {

    /**
     * A struct namespace for methods to handle endianness in integer values.
     */
    template <class IntType>
    struct endianness {
    public:
        
        /// Converts from an integer in the native representation in the machine
        /// architecture to a big-endian representation
        static inline IntType to_big_endian(IntType value);
        
        /// Converts from a big-endian integer to the native representation in
        /// the machine architecture
        static inline IntType from_big_endian(IntType value);
        
    private:
        
        /// Returns the integer in the opposite endianness representation it currently
        /// has
        static inline IntType swap_endianness(IntType value);
        
        /// Returns true if the architecture is big-endian, otherwise false
        static inline bool arch_is_big_endian();
        
        static_assert(std::is_integral<IntType>::value,
                      "Endianness only applies to integer data types");
    };
    
    ////////////////////////////
    /// Template implementations
    ////////////////////////////
    
    
    template <class IntType>
    inline IntType endianness<IntType>::to_big_endian(IntType value) {
        return arch_is_big_endian() ? value : swap_endianness(value);
    }
    
    template <class IntType>
    inline IntType endianness<IntType>::from_big_endian(IntType value) {
        // these turn out to be identical functions, but having both aliases still
        // seems cognitively useful
        return to_big_endian(value);
    }
    
    template <class IntType>
    inline IntType endianness<IntType>::swap_endianness(IntType value) {
        
        IntType swapped;
        
        uint8_t* from = (uint8_t*) &value;
        uint8_t* to = (uint8_t*) &swapped;
        
        for (int i = 0; i < sizeof(IntType); ++i) {
            to[i] = from[sizeof(IntType) - i - 1];
        }
        
        return swapped;
    }
    
    // TODO: this method will not detect endianness correctly on 1-byte
    // integers, but endianness is irrelevant for them anyway...
    template <class IntType>
    inline bool endianness<IntType>::arch_is_big_endian() {
        
        // mark volatile so the compiler won't optimize it away
        volatile IntType val = 1;
        
        uint8_t* bytes = (uint8_t*) &val;
        
        // the 1 is only set at the lowest memory address if the architecture
        // is little-endian
        return !bytes[0];
    }
}

#endif
