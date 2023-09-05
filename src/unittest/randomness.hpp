#ifndef VG_UNITTEST_RANDOMNESS_HPP_INCLUDED
#define VG_UNITTEST_RANDOMNESS_HPP_INCLUDED

/**
 * \file randomness.hpp
 * Provide C++11 <random>-compatible randomness that is seedable via Catch.hpp RNG seeding.
 *
 * NO TEST SHOULD USE RANDOMNESS FROM ANY OTHER SOURCE!
 */

#include <random>
#include <cstdlib>

namespace vg {
namespace unittest {

/**
 * Return a random seed for testing purposes.
 * Respects catch.hpp --rng-seed option to `vg test`.
 *
 * USE INSTEAD OF A NEW std::random_device FOR ALL TESTING PURPOSES!
 */
inline unsigned int test_seed_source() {
    // TODO: make thread safe
    return (unsigned int) rand();
}

}
}

#endif
