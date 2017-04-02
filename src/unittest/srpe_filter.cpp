/** \file
 * unittest/sift_test.cpp: tests for both paired and unpaired read filtering
 */

#include "catch.hpp"
#include "filter.hpp"

namespace vg
{
namespace unittest
{

using namespace std;

TEST_CASE("paired-end reads that are one-end anchored are correctly reported.", "[scrub]")
{

    // Make some deletion reads

    // make some insertion reads

    // make some inversion reads
    Filter ff;
    const string oea_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {"name": "x", "mapping" : []},
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";

    const string oea_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {},
                        "mapping_quality" : 0,
                        "fragment_next" : "mate_one"
                }
                )";

    Alignment uno;
    Alignment dos;
    json2pb(uno, oea_one_json.c_str(), oea_one_json.size());
    json2pb(dos, oea_two_json.c_str(), oea_two_json.size());

    const string norm_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {"name": "x", "mapping" : []},
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
    const string norm_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {},
                        "mapping_quality" : 1000,
                        "fragment_next" : "mate_one"
                }
                )";

    Alignment uno_normal;
    Alignment dos_normal;
    json2pb(uno_normal, norm_one_json.c_str(), norm_one_json.size());
    json2pb(dos_normal, norm_two_json.c_str(), norm_two_json.size());
    pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno_normal, dos_normal);

    SECTION("Read pairs that are one-end-anchored are kept in standard mode.")
    {

        pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno, dos);

        // It has to say it didn't
        // trim it
        REQUIRE(ret.first.name() == "mate_one");
        REQUIRE(ret.second.name() == "mate_two");
        // And it has to not mess it up
    }

    SECTION("Empty alignments are returned for reads that are not one-end-anchored in standard mode.")
    {

        // It has to say it didn't
        // trim it
        REQUIRE(ret.first.name() == "");
        REQUIRE(ret.second.name() == "");
        // And it has to not mess it up
    }
    SECTION("Read pairs that are one-end-anchored are removed in inverse mode.")
    {
        ff.set_inverse(true);
        pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno, dos);

        // It has to say it didn't
        // trim it
        REQUIRE(ret.first.name() == "");
        REQUIRE(ret.second.name() == "");
        // And it has to not mess it up
    }

    SECTION("Read pairs that are not one-end-anchored are returned in inverse mode.")
    {
        ff.set_inverse(true);
        pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno_normal, dos_normal);

        // It has to say it didn't
        // trim it
        REQUIRE(ret.first.name() == "mate_one");
        REQUIRE(ret.second.name() == "mate_two");
        // And it has to not mess it up
    }
}
}
}
