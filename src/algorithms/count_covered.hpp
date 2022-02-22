#ifndef VG_ALGORITHMS_COUNT_COVERED_HPP_INCLUDED
#define VG_ALGORITHMS_COUNT_COVERED_HPP_INCLUDED

/**
 * \file count_covered.hpp
 *
 * Sweep-line algorithm to count the number of positions covered by a set of
 * intervals.
 */

namespace vg {
namespace algorithms {

using namespace std;

/**
 * Count, from begin to end, the number of positions covered by ranges in the
 * given collection. The collection will be sorted in place.
 *
 * The collection must be a have a begin(), end(), and random [] access (like a
 * vector). 
 *
 * No boundaries are needed because no positions can be covered without
 * segments representing them.
 */
template<typename Collection>
size_t count_covered(Collection& segments) {

    if (segments.empty()) {
        // Protect against no segments
        return 0;
    }

    std::sort(segments.begin(), segments.end());
    auto curr_begin = segments[0].first;
    auto curr_end = segments[0].second;
    
    size_t total = 0;
    for (size_t i = 1; i < segments.size(); i++) {
        if (segments[i].first >= curr_end) {
            total += (curr_end - curr_begin);
            curr_begin = segments[i].first;
            curr_end = segments[i].second;
        }
        else if (segments[i].second > curr_end) {
            curr_end = segments[i].second;
        }
    }
    total += (curr_end - curr_begin);
    return total;
}

}
}

#endif

