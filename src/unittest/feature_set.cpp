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
    
    // Run the BED through the load
    FeatureSet features;
    features.load_bed(in);
    
    // Look at the features that are in there
    auto& loaded = features.get_features("seq1");
    REQUIRE(loaded.size() == 1);
    REQUIRE(loaded[0].path_name == "seq1");
    REQUIRE(loaded[0].first == 5);
    REQUIRE(loaded[0].last == 10);
    REQUIRE(loaded[0].feature_name == "record");
    REQUIRE(loaded[0].extra_data.size() == 0);
    
    // Run through the save
    features.save_bed(out);
    
    // Make sure it didn't change anything
    REQUIRE(out.str() == "seq1\t5\t10\trecord\n");

}

TEST_CASE("Internal deletions move the end left", "[featureset][simplify]") {

    // Make a BED stream to read
    stringstream in("seq1\t5\t10\trecord\n");
    // And a BED stream to write to
    stringstream out;
    
    // Run the BED through the load
    FeatureSet features;
    features.load_bed(in);
    
    // Add an internal deletion
    features.on_path_edit("seq1", 6, 1, 0);
    
    // Run through the save
    features.save_bed(out);
    
    // Make sure it moved the end left
    REQUIRE(out.str() == "seq1\t5\t9\trecord\n");

}

}
}

