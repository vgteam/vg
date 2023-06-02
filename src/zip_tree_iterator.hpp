#ifndef VG_ZIP_TREE_ITERATOR_HPP_INCLUDED
#define VG_ZIP_TREE_ITERATOR_HPP_INCLUDED

#include <vector>
#include <utility>
#include <cstdint>
#include "position.hpp"

/** \file zip_tree_iterator.hpp
 * Iterator for querying predecessors of zipcoded positions in a graph space.
 */

namespace vg{
using namespace std;



/**
 * Interface for something that represents a tree of positions and their
 * distances, based on their zipcodes. Tree would be constructed from pairs of
 * positions and zipcodes, and stored in a mostly-linear representation as
 * specified in
 * https://github.com/benedictpaten/long_read_giraffe_chainer_prototype/blob/b590c34055474b0c901a681a1aa99f1651abb6a4/zip_tree_iterator.py.
 */
class ZipTreeStringInterface {
public:

    class enum Bound {
        SNARL_START,
        SNARL_END,
        CHAIN_START,
        CHAIN_END,
    };

    union Value {
        pos_t position;
        size_t distance;
        ZipTreeBound bound;
    };

    class enum ValueType {
        POSITION,
        DISTANCE,
        BOUND
    };

    /**
     * Base class for a simple iterator that looks left one step in the string.
     */
    class ReverseIterator {
        virtual ~ZipTreeStringIterator() = default;

        /// Move one left
        virtual ZipTreeStringIterator& operator++();

        /// Compare for equality to see if we hit end (the past-the-left position)
        virtual bool operator==(const ZipTreeStringIterator& other) const;

        /// Compare for inequality
        inline bool operator!=(const ZipTreeStringIterator& other) const {
            return !(*this == other);
        }

        /// Produce a type-tagged union expressing the value that the iterator is at in the string.
        std::pair<ValueType, Value>* operator*() const;
    };

};

    
    
}
#endif
