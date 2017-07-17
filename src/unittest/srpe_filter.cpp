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
    Alignment clean_first;
    Alignment clean_fwd_first;
    Alignment second;
    Alignment second_rev;
    Alignment f_rev;
    Alignment s_forward;
    Alignment s_flip;
    Alignment f_flip;
    json2pb(clean_first, one_json.c_str(), one_json.size());
    json2pb(clean_fwd_first, one_json.c_str(), one_json.size());
    json2pb(second, two_json.c_str(), two_json.size());
    json2pb(second_rev, two_json.c_str(), two_json.size());
    json2pb(f_rev, one_rev_json.c_str(), one_rev_json.size());
    json2pb(s_forward, two_fwd_json.c_str(), two_fwd_json.size());
    json2pb(s_flip, two_flip_json.c_str(), two_flip_json.size());
    json2pb(f_flip, one_flip_json.c_str(), one_flip_json.size());

    SECTION("Paired_orientation_filter ignores normal pairs and catches flipped ones"){
        bool norm = ff.pair_orientation_filter(clean_first, second);
        bool rev_rev = ff.pair_orientation_filter(f_rev, second_rev);
        bool fwd_fwd = ff.pair_orientation_filter(clean_fwd_first, s_forward);
        bool flippies = ff.pair_orientation_filter(f_flip, s_flip);
     
        REQUIRE(norm == false);
        REQUIRE(rev_rev);
        REQUIRE(fwd_fwd);
        REQUIRE(flippies);

        REQUIRE(clean_first.name() != "");
        REQUIRE(clean_first.read_on_reverse_strand() == false);
        REQUIRE(clean_first.mate_on_reverse_strand() == true);
        REQUIRE(second.read_on_reverse_strand() == true);
        REQUIRE(second.mate_on_reverse_strand() == false);

        REQUIRE(f_rev.name() != "");
        REQUIRE(f_rev.read_on_reverse_strand() == true);
        REQUIRE(f_rev.mate_on_reverse_strand() == true);
        REQUIRE(second_rev.read_on_reverse_strand() == true);
        REQUIRE(second_rev.mate_on_reverse_strand() == true);

        REQUIRE(clean_fwd_first.name() != "");
        REQUIRE(clean_fwd_first.name() != "");
        REQUIRE(clean_fwd_first.read_on_reverse_strand() == false);
        REQUIRE(clean_fwd_first.mate_on_reverse_strand() == false);
        REQUIRE(s_forward.read_on_reverse_strand() == false);
        REQUIRE(s_forward.mate_on_reverse_strand() == false);

        REQUIRE(f_flip.name() != "");
        REQUIRE(s_flip.name() != "");
        REQUIRE(f_flip.read_on_reverse_strand() == false);
        REQUIRE(f_flip.mate_on_reverse_strand() == true);
        REQUIRE(s_flip.read_on_reverse_strand() ==  true);
        REQUIRE(s_flip.mate_on_reverse_strand() == false);
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
        bool t = ff.soft_clip_filter(first);
        REQUIRE(t == true);
        REQUIRE(first.soft_clipped() == true);
    }
    SECTION("Soft clip filter ignores non-clipped reads"){
        REQUIRE(ff.soft_clip_filter(second) == false);
        ff.set_soft_clip_limit(1000);
        REQUIRE(!ff.soft_clip_filter(first));
        REQUIRE(first.soft_clipped() == false);
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


    SECTION("Soft clip trimmer removes clipped portions"){
        Alignment a = ff.remove_clipped_portion(first);
        REQUIRE(a.soft_clipped() == false);
        // Should have no side effects on the input alignment, even though our alignment is indeed softclipped
        REQUIRE(first.soft_clipped() == false);
        REQUIRE(a.name() != "");
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
        REQUIRE(ff.unmapped_filter(first) == true);

    }
    SECTION("Unmapped filter ignores reads with no score"){
            REQUIRE(ff.unmapped_filter(second) == false);

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
    bool oeas = ff.one_end_anchored_filter(first, second);
    REQUIRE(oeas == true);
    REQUIRE(first.read_mapped() == false);
    REQUIRE(first.mate_unmapped() == false);

}

}
}
