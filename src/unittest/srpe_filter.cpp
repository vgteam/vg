/** \file
 * unittest/sift_test.cpp: tests for both paired and unpaired read filtering
 */

#include "catch.hpp"
#include "filter.hpp"
#include "../utility.hpp"
#include "../json2pb.h"

namespace vg
{
namespace unittest
{

TEST_CASE("filter correctly handles clipped/unclipped reads.", "[sift]") {

    Filter ff;
    ff.set_soft_clip_limit(10);

    const string one_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCA"}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    const string two_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    Alignment first;
    Alignment second;
    json2pb(first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());

    SECTION("Soft clip filter catches clipped reads"){
    REQUIRE(ff.soft_clip_filter(first).name() != "");
    }
    SECTION("Soft clip filter ignores non-clipped reads"){
        REQUIRE(ff.soft_clip_filter(second).name() == "");
        ff.set_soft_clip_limit(1000);
        REQUIRE(ff.soft_clip_filter(first).name() == "");
    }
}

TEST_CASE("Softclipped portions can be trimmed from reads.", "[sift]") {

    Filter ff;
    ff.set_soft_clip_limit(10);

    const string one_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCA"}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    const string two_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    Alignment first;
    Alignment second;
    json2pb(first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());

    SECTION("Soft clip filter removes clipped portions"){
        Alignment a = ff.remove_clipped_portion(first);
        cout << pb2json(a);
        REQUIRE(ff.soft_clip_filter(a).name() == "");
    }

}

TEST_CASE("Unmapped reads are caught by filter and mapped ones ignored.", "[sift]"){

    Filter ff;

    const string one_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                  "name": "Random-147892/2",
                "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                "fragment_prev": {"name": "Random-147892/1"}}
                )";

    const string two_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60,
                     "score": 79,
                      "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    Alignment first;
    Alignment second;
    json2pb(first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());


    SECTION("Unmapped filter catches reads with no score"){
    REQUIRE(ff.unmapped_filter(first).name() != "");

    }
    SECTION("Unmapped filter ignores reads with no score"){
            REQUIRE(ff.unmapped_filter(second).name() == "");

    }


}

}
}
