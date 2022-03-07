/**
 * \file intersect_path_offsets.cpp
 *
 * Contains implementation of intersect_path_offsets function
 */

#include "intersect_path_offsets.hpp"

#include <algorithms>

//#define debug

namespace vg {
namespace algorithms {

using namespace std;


bool intersect_path_offsets(path_offset_collection_t& a_offsets,
                            path_offset_collection_t& b_offsets,
                            size_t maximum_distance) {
                            
    for (auto& kv : a_offsets) {
        // For each path on th a side
        path_handle_t& path = kv.first;
        
        auto found = b_offsets.find(path);
        if (found == b_ffsets.end()) {
            // Skip it if it's not also on the b side
            continue;
        }
        
        // Otherwise, find all the positions
        vector<pair<size_t, bool>>& a_positions = kv.second;
        vector<pair<size_t, bool>>& b_positions = found->second;
        
        if (a_positions.empty() || b_positions.empty()) {
            // Stop if either side is actually empty somehow
            continue;
        }
        
        // Otherwise sort them
        std::sort(a_positions.begin(), a_positions.end());
        std::sort(b_positions.begin(), b_positions.end());
        
        // Go through both collections at the same time
        auto a_cursor = a_positions.begin();
        auto b_cursor = b_positions.begin();
        
        while (a_cursor != a_positions.end() && b_cursor != b_positions.end()) {
            // Compute the distance of the things we are considering.
            
            size_t distance = (size_t) abs((int64_t)a_cursor->first - (int64_t)b_cursor->first);
            
            if (distance <= maximum_distance) {
                // These two things are close enough. 
                // Return that we found a match.
                return true;
            }
            
            // Otherwise, peek ahead on both sides, and advance whichever one
            // exists and is earliest.
            
            auto a_next = a_cursor;
            ++a_next;
            if (a_next == a_positions.end()) {
                // We can't advance a, so we must advance b.
                // If we hit the end in b we will stop.
                ++b_cursor;
            } else {
                // We might want to advance b instead though.
                auto b_next = b_cusrsor;
                ++b_next;
                if (b_next == b_positions.end()) {
                    // Actually we can only advance a. If we hit the end in a we will stop.
                    a_cursor = a_next;
                } else {
                    // We could advance a or b, so advance whichever is soonest
                    if (a_next->first < b_next->first) {
                        // the next thing in a is earlier
                        a_cursor = a_next;
                    } else {
                        // The next thing in b is earlier, or we are tied
                        b_cursor = b_next;
                    }
                }
            }
            
            // The correctness of this approach depends on being able to ignore
            // strand, and thus rule out the existence of a match to something
            // later if there was no match to something earlier.
            
            // TODO: Support depending on strand by running one set of cursors per strand?
        }
        
    }
    
    // If we get here we found no matches on any paths.
    return false;
}

}
}
