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

TEST_CASE("GAMindex binning works on a large number", "[gam][gamindex]") {

    for (id_t to_bin : {0ULL, 1ULL, 10ULL, 0xFFFFFFFFFFFFFFFFULL, 0xFACEDEADCAFEBEEFULL}) {

        auto bins = GAMIndex::bins_of_id(to_bin);

        // We need one bin per bit for all bits but the last one, plus one top-level 0 bin for everything.
        REQUIRE(bins.size() == CHAR_BIT * sizeof(id_t));

        // The bins should end with the least specific bin (0)
        REQUIRE(bins.back() == 0);
        
        // Bin levels go from 0 to bits-1.
        // The first bin should be 2^(bits - 1) - 1 + id >> 1
        REQUIRE(bins.front() == ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1 + ((GAMIndex::bin_t)to_bin >> 1));
        
        // The first bin is the most specific bin, which is the bin of the number and itself.
        REQUIRE(bins.front() == GAMIndex::common_bin(to_bin, to_bin));
    }
    
}

TEST_CASE("GAMindex binning works on adjacent numbers", "[gam][gamindex]") {
    
    // The common bin of an even and the next odd number should be the two numbers right shifted by one, pulss an offset.
    // For an odd and the next even, it should be them shifted by 2, pluss a smaller offset.
    
    // what offsets do we expect for bins based on their size?
    auto SIZE_2_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 1)) - 1;
    auto SIZE_4_OFFSET = ((GAMIndex::bin_t)0x1 << (CHAR_BIT * sizeof(id_t) - 2)) - 1;
    
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

TEST_CASE("GAMindex can look up inserted ranges", "[gam][gamindex]") {
    // Make an empty index
    GAMIndex index;

    // Add some ID-sorted groups
    index.add_group(1, 5, 0, 100);
    index.add_group(3, 7, 100, 200);
    index.add_group(6, 9, 200, 300);
    // Being sorted by lowest ID doesn't mean you are always sorted by highest ID
    index.add_group(7, 8, 300, 400);
    index.add_group(100, 110, 400, 500);
    
    // Look for node 1
    auto found = index.find(1);
    // We should find the 1 run from 0 to 100
    REQUIRE(found.size() == 1);
    REQUIRE(found[0].first == 0);
    REQUIRE(found[0].second == 100);
    
    // Look for node 7
    found = index.find(7);
    // We should find one combined run from 100 to 400
    REQUIRE(found.size() == 1);
    REQUIRE(found[0].first == 100);
    REQUIRE(found[0].second == 400);
    
    // Look for node 999 which nothing can touch
    found = index.find(999);
    REQUIRE(found.size() == 0);
    
    
}


}
}
