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

TEST_CASE("filter correctly handles inverted read pairs", "[sift"){
    Filter ff;
    const string one_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": false}, "edit": [{"to_length": 46, "sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCA"}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
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

    const string one_rev_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCA"}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";
    const string two_fwd_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": false}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";

    const string one_flip_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": false}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": -1000}]}
                )";
    const string two_flip_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60, "score": 79, "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": -1000}]}
                )";
    Alignment first;
    Alignment second;
    Alignment f_rev;
    Alignment s_forward;
    Alignment s_flip;
    Alignment f_flip;
    json2pb(first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());
    json2pb(f_rev, one_rev_json.c_str(), one_rev_json.size());
    json2pb(s_forward, two_fwd_json.c_str(), two_fwd_json.size());
    json2pb(s_flip, two_flip_json.c_str(), two_flip_json.size());
    json2pb(f_flip, one_flip_json.c_str(), one_flip_json.size());

    SECTION("Paired_orientation_filter ignores normal pairs"){
        pair<Alignment, Alignment> norm = ff.pair_orientation_filter(first, second);
        pair<Alignment, Alignment> rev_rev = ff.pair_orientation_filter(f_rev, second);
        pair<Alignment, Alignment> fwd_fwd = ff.pair_orientation_filter(first, s_forward);
        pair<Alignment, Alignment> flippies = ff.pair_orientation_filter(f_flip, s_flip);
     
        REQUIRE(norm.first.name() == "");
        REQUIRE(norm.second.name() == "");

        REQUIRE(rev_rev.first.name() != "");
        REQUIRE(rev_rev.first.name() != "");
        REQUIRE(rev_rev.first.read_on_reverse_strand() == true);
        REQUIRE(rev_rev.first.mate_on_reverse_strand() == true);
        REQUIRE(rev_rev.second.read_on_reverse_strand() == true);
        REQUIRE(rev_rev.second.mate_on_reverse_strand() == true);

        REQUIRE(fwd_fwd.first.name() != "");
        REQUIRE(fwd_fwd.first.name() != "");
        REQUIRE(fwd_fwd.first.read_on_reverse_strand() == false);
        REQUIRE(fwd_fwd.first.mate_on_reverse_strand() == false);
        REQUIRE(fwd_fwd.second.read_on_reverse_strand() == false);
        REQUIRE(fwd_fwd.second.mate_on_reverse_strand() == false);

        REQUIRE(flippies.first.name() != "");
        REQUIRE(flippies.first.name() != "");
        REQUIRE(flippies.first.read_on_reverse_strand() == false);
        REQUIRE(flippies.first.mate_on_reverse_strand() == true);
        REQUIRE(flippies.second.read_on_reverse_strand() ==  true);
        REQUIRE(flippies.second.mate_on_reverse_strand() == false);
    }
}


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
        REQUIRE(ff.soft_clip_filter(first).soft_clipped() == true);
    }
    SECTION("Soft clip filter ignores non-clipped reads"){
        REQUIRE(ff.soft_clip_filter(second).name() == "");
        ff.set_soft_clip_limit(1000);
        REQUIRE(ff.soft_clip_filter(first).name() == "");
        REQUIRE(ff.soft_clip_filter(first).soft_clipped() == false);
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
        REQUIRE(ff.soft_clip_filter(a).name() == "");
        REQUIRE(ff.soft_clip_filter(a).soft_clipped() == false);
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

TEST_CASE("One-end-anchored filter works", "[sift]"){
    const string one_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                  "name": "Random-147892/2",
                "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                "fragment_prev": {"name": "Random-147892/1"}
                }
                )";

    const string two_json = R"(
                {"sequence": "TTTCGCAGCGATCGCGACCTTAGCAAAGAGTCCCCTGCTAGTTCCAGGCATTTGGCTGCTGCTCAGAATCTACGGGCTATATTAAAGCAACGGGAGACTTGAGTCCAATCCGTTACCTGTAGTCA",
                 "path": {"mapping": [{"position": {"node_id": 120, "offset": 866, "is_reverse": true}, "edit": [{"to_length": 46, "from_length": 46}, {"from_length": 79, "to_length": 79}], "rank": 1}]},
                  "name": "Random-147892/2",
                   "quality": "ISIiIiImJiYeJiYfJiYmJhAcISYmJiYmJiYmJiYmJiYmJiYmJiYmJiYmJg4mJiYmJiYmECYmJiYmJiImIiYmJiUmJiYmJiYmHSYmJiQmJiYmJiYmJiYlJSUmJiYmJSYmJiYiJiYmJiYmJiYlJhwmJiYmJiYPJiMmJiYmJiY=",
                    "mapping_quality": 60,
                     "score": 60,
                      "fragment_prev": {"name": "Random-147892/1"},
                     "identity": 0.63200000000000001, "fragment": [{"name": "Random", "length": 1000}]}
                )";
            Alignment first;
    Alignment second;
    json2pb(first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());

    Filter ff;
    pair<Alignment, Alignment> oeas = ff.one_end_anchored_filter(first, second);
    REQUIRE(oeas.first.name() != "");
    REQUIRE(oeas.first.read_mapped() == false);
    REQUIRE(oeas.first.mate_unmapped() == false);

}

}
}
