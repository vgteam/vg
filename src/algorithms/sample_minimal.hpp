#ifndef VG_ALGORITHMS_SAMPLE_MINIMAL_HPP_INCLUDED
#define VG_ALGORITHMS_SAMPLE_MINIMAL_HPP_INCLUDED

/**
 * \file
 * Minimizer (sub)sampling algorithm, as explained in the Winnowmap paper, Jain et al. 2020.
 * Goes through read space and samples all candidates that are minimal in a sliding window of a given size.
 */
 
#include <functional>

namespace vg {
namespace algorithms {

using namespace std;


/**
 * Sample the minimal elements in windows of the given size. Uses get_bounds to
 * get inclusive-start, exclusive-end coordinates for elements. Uses
 * should_beat to compare elements. If an element is minimal for a window,
 * calls sample for that element.
 *
 * You can use should_beat to control tie behavior. If it acts as a less-than
 * comparator, and returns false for ties, tied elements will all be sampled.
 * If it acts as less-than-or-equal-to, and returns true for ties, the
 * latest-occurring element will be sampled in case of ties.
 *
 * Elements must be sorted by start and all the same length.
 *
 * Unlike the minimizer sampling algorithm given in Jain et al. 2020., we have
 * to make sure to support multiple elements on the same start position, and
 * zero elements on some start positions.
 *
 * sample will be called at least once for each element minimal in some window.
 * It will not necessarily be called once per window.
 */
void sample_minimal(size_t count, size_t element_length, size_t window_size, size_t sequence_length, const std::function<size_t(size_t)>& get_start, const std::function<bool(size_t, size_t)>& should_beat, const std::function<void(size_t)>& sample);

}

}

#endif
