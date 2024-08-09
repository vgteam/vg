#ifndef VG_ALGORITHMS_PAD_BAND_HPP_INCLUDED
#define VG_ALGORITHMS_PAD_BAND_HPP_INCLUDED

/**
 * \file pad_band.hpp
 *
 * Defines algorithm for computing band padding for banded alignment.
 */

#include "../handle.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace algorithms {

using namespace std;

/// Get a band padding function that scales with the expected distance of a random
/// walk, memoized out to the given length.
std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_random_walk(double band_padding_multiplier = 1.0, size_t band_padding_memo_size = 2000);

/// Get a band padding function that scales the expected distance of a random
/// walk, memoized out to the given length, using the minimum of graph size and
/// read size as the length.
std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_min_random_walk(double band_padding_multiplier = 1.0, size_t band_padding_memo_size = 2000);

/// Get a band padding function that uses a constant value.
std::function<size_t(const Alignment&, const HandleGraph&)> pad_band_constant(size_t band_padding);

}
}

#endif
