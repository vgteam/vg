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

TEST_CASE("paired-end reads that are one-end anchored are correctly reported.", "[sift]")
{

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

  

}
}
}
