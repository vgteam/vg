///
///  \file genome_state.cpp
///
///  Unit tests for the GAMIndex which indexes seekable GAM files by node ID
///

#include <iostream>
#include "catch.hpp"
#include "../gam_index.hpp"


namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("GAMindex windowing works correctly", "[gam][gamindex]") {

    for (size_t i = 0; i < 10000000; i += 83373) {
        REQUIRE(GAMIndex::window_of_id(i) == i / 256);
    }

}

TEST_CASE("GAMindex binning works correctly", "[gam][gamindex]") {

    auto bins = GAMIndex::bins_of_id(0xACAC);

    // We need one bin per bit for all bits but the last one.
    REQUIRE(bins.size() == CHAR_BIT * sizeof(id_t) - 1);

    // The bins should end with the least specific bin (0)
    REQUIRE(bins.back() == 0);
    
    // what offsets do we expect for bins based on their size?
    auto SIZE_2_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1;
    auto SIZE_4_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 2)) - 1;
    
    // The first bin should be 2^(bits - 1) - 1 + id >> 1
    REQUIRE(bins.front() == ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1 + ((GAMIndex::bin_t)0xACAC >> 1));
    
    // The common bin of an even and the next odd number should be the two numbers right shifted by one, pulss an offset.
    // For an odd and the next even, it should be them shifted by 2, pluss a smaller offset.
    
    for (size_t i = -10; i < 10000; i++) {
        auto bin_found = GAMIndex::common_bin(i, i + 1);
        auto bin_found2 = GAMIndex::common_bin(i + 1, i);
        
        // Should work in any order
        REQUIRE(bin_found == bin_found2);
        
        if (i == -1) {
            // The common bin between -1 and 0 has to be bin 0, because that's where the discontinuity falls.
            // Not a problem because we don't use negative node IDs in real life.
            REQUIRE(bin_found == 0);
        } else {
            if (i % 2 == 0) {
                // Even number and next odd
                REQUIRE(bin_found == SIZE_2_OFFSET + ((GAMIndex::bin_t)i >> 1));
            } else {
                // Odd number and next even
                REQUIRE(bin_found == SIZE_4_OFFSET + ((GAMIndex::bin_t)i >> 2));
            }
        }
    }
    

}


}
}
