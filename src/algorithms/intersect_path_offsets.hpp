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
 * determine if any positions in the two sets are on the
 * same path, within the given maximum distance.
 *
 * The set expected to visit fewer paths should be passed first.
 *
 * Orientation is ignored.
 *
 * See nearest_offsets_in_paths() for a function that generates these.
 *
 * Needs write access to its input position collections, so it can sort them in
 * place.
 */
bool intersect_path_offsets(path_offset_collection_t& a_offsets,
                            path_offset_collection_t& b_offsets,
                            size_t maximum_distance);
    

}
}

#endif

