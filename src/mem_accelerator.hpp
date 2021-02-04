/**
 * \file mem_accelerator.hpp
 *
 * Defines an index for accelerating GCSA2 queries
 */

#ifndef VG_MEM_ACCELERATOR_HPP_INCLUDED
#define VG_MEM_ACCELERATOR_HPP_INCLUDED

#include <cstdint>
#include <string>
#include <gcsa/gcsa.h>
#include <sdsl/int_vector.hpp>

namespace vg {

using namespace std;

/*
 * An auxilliary index that accelerates the initial steps of
 * MEM-finding in
 */
class MEMAccelerator {
public:
    
    MEMAccelerator() = default;
    MEMAccelerator(const gcsa::GCSA& gcsa_index, size_t k);
    
    // return the length of k-mers that are memoized
    inline int64_t length() const;
    
    // look up the GCSA range that corresponds to a k-length
    // string ending at the indicated position. client code
    // is responsible for ensuring that the string being
    // accessed is at least length k and consists only of ACGT
    // characters
    gcsa::range_type memoized_LF(string::const_iterator last) const;
    
private:
    
    inline int64_t encode(char c) const;
    
    // the size k-mer we'll index
    const int64_t k = 1;
    // the actual table
    sdsl::int_vector<> range_table;
    
};

inline int64_t MEMAccelerator::length() const {
    return k;
}

inline int64_t MEMAccelerator::encode(char c) const {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

}

#endif
