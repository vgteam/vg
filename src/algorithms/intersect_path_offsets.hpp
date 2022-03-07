#ifndef VG_ALGORITHMS_INTERSECT_PATH_OFFSETS_HPP_INCLUDED
#define VG_ALGORITHMS_INTERSECT_PATH_OFFSETS_HPP_INCLUDED

/**
 * \file intersect_path_offsets.hpp
 *
 * Defines algorithm for finding whether any offsets on paths in one set are
 * near any offsets on paths in a different set.
 */

#include "../handle.hpp"

#include "nearest_offsets_in_paths.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    
/**
 * Given two maps from path handle to (position, orientation) pair vectors,
 * determine if any positions in the two sets are on the same path, within the
 * given maximum distance.
 *
 * The set expected to have more visits should be passed first.
 *
 * Orientation is ignored.
 *
 * The first set must be sorted, for binary search. We run binary search for
 * each item in the second set, so the first set should be the larger one.
 *
 * We run in b log a time.
 */
bool intersect_path_offsets(const path_offset_collection_t& a_offsets,
                            const path_offset_collection_t& b_offsets,
                            size_t maximum_distance);
                            
/**
 * Sort path offsets, so intersect_path_offsets() can use them as a target.
 */
void sort_path_offsets(path_offset_collection_t& offsets);
    

}
}

#endif

