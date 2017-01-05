/** \file
 *
 * Unit tests for the FeatureSet, which keeps some BED-like regions in sync with
 * changes to a graph.
 */

#include <iostream>
#include <sstream>
#include "../feature_set.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("FeatureSet can load and save a BED record", "[featureset][simplify]") {

    // Make a BED stream to read
    stringstream in("seq1\t5\t10\trecord\n");
    // And a BED stream to write to
    stringstream out;
    
    // Run the BED through the load and save
    FeatureSet features;
    features.load_bed(in);
    features.save_bed(out);
    
    // Make sure it didn't change anything
    REQUIRE(out.str() == "seq1\t5\t10\trecord\n");

}

}
}

