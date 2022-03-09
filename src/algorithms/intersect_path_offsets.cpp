/**
 * \file intersect_path_offsets.cpp
 *
 * Contains implementation of intersect_path_offsets function
 */

#include "intersect_path_offsets.hpp"

#include <algorithm>

//#define debug

namespace vg {
namespace algorithms {

using namespace std;


bool intersect_path_offsets(const path_offset_collection_t& a_offsets,
                            const path_offset_collection_t& b_offsets,
                            size_t maximum_distance) {
                            
    for (auto& kv : b_offsets) {
        // For each path on the b side
        const path_handle_t& path = kv.first;
        
        auto found = a_offsets.find(path);
        if (found == a_offsets.end()) {
            // Skip it if it's not also on the a side
            continue;
        }
        // If it is on the a side, we'll do search against this sorted list.
        auto& target_positions = found->second;
        
        for (auto& b_position : kv.second) {
            // For each offset on the b side, we need to do a binary search.
            
            // Find the nearest thing to our right, if any.
            auto falling_after_it = std::lower_bound(target_positions.begin(), target_positions.end(), b_position);
            if (falling_after_it != target_positions.begin()) {
                // There's also going to be something to our left. Check that first.
                auto falling_before_it = falling_after_it;
                --falling_before_it;
                if (b_position.first - falling_before_it->first <= maximum_distance) {
                    // It is close enough before to constitute a hit
                    return true;
                }
            }
            if (falling_after_it != target_positions.end()) {
                // The thing to our right actually exists. Check it too.
                if (falling_after_it->first - b_position.first <= maximum_distance) {
                    // It is close enough to constitute a hit
                    return true;
                }
            }
        }
    }
        
    // If we get here we found no matches on any paths.
    return false;
}

void sort_path_offsets(path_offset_collection_t& offsets) {
    for (auto& kv : offsets) {
        // Sort along each path.
        std::sort(kv.second.begin(), kv.second.end());
    }
}

}
}
