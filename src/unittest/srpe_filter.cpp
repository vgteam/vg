/**
 * unittest/readfilter.cpp: test cases for the read filtering/transforming logic
 */

#include "catch.hpp"
#include "filter.hpp"

namespace vg {
    namespace unittest {

        using namespace std;

        TEST_CASE("paired-end reads that are one-end anchored are correctly reported.", "[scrub]") {

            Filter ff;


            SECTION("Read pairs that are one-end-anchored are kept in standard mode.") {
                const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {"name": "x", "mapping" : []},
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
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
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "mate_one");
                REQUIRE(ret.second.name() == "mate_two");
                // And it has to not mess it up


            }

            SECTION("Empty alignments are returned for reads that are not one-end-anchored in standard mode.") {
                const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {"name": "x", "mapping" : []},
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {},
                        "mapping_quality" : 1000,
                        "fragment_next" : "mate_one"
                }
                )";

                Alignment uno;
                Alignment dos;
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "");
                REQUIRE(ret.second.name() == "");
                // And it has to not mess it up


            }
            SECTION("Read pairs that are one-end-anchored are removed in inverse mode.") {
                ff.set_inverse(true);

                const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {"name": "x", "mapping" : []},
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
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
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.one_end_anchored_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "");
                REQUIRE(ret.second.name() == "");
                // And it has to not mess it up


            }

            SECTION("Discordant orientation filter does not flag reads aligned normally."){
                ff.set_inverse(false);
                const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                        {"position": {"node_id": 0, "is_reverse" : true}}
                        ]
                        },
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "mapping_quality" : 100,
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                {"position": {"node_id": 0, "is_reverse" : false}}
                            ]
                        },
                        "fragment_next" : "mate_one"
                }
                )";

                Alignment uno;
                Alignment dos;
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.orientation_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "");
                REQUIRE(ret.second.name() == "");
                // And it has to not mess it up



            }

            SECTION("Discordant orientation filter flags reads that look like --> -->."){
                const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                        {"position": {"node_id": 0, "is_reverse" : false}}
                        ]
                        },
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "mapping_quality" : 100,
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                {"position": {"node_id": 0, "is_reverse" : false}}
                            ]
                        },
                        "fragment_next" : "mate_one"
                }
                )";

                Alignment uno;
                Alignment dos;
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.orientation_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "mate_one");
                REQUIRE(ret.second.name() == "mate_two");
            }

            SECTION("Discordant orientation filter flags reads that look like <-- <--"){
            const string mate_one_json = R"(
                {
                    "name" : "mate_one",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                        {"position": {"node_id": 0, "is_reverse" : true}}
                        ]
                        },
                        "fragment_next" : "mate_two",
                        "mapping_quality": 100
                }
                )";
                const string mate_two_json = R"(
                {
                    "name" : "mate_two",
                        "sequence" : "AATTGGCCAAGGCA",
                        "quality" : "MDAwMDAwMDAwMDAwMDAw",
                        "mapping_quality" : 100,
                        "path" : {
                            "name": "x",
                            "mapping" : [
                                {"position": {"node_id": 0, "is_reverse" : true}}
                            ]
                        },
                        "fragment_next" : "mate_one"
                }
                )";

                Alignment uno;
                Alignment dos;
                json2pb(uno, mate_one_json.c_str(), mate_one_json.size());
                json2pb(dos, mate_two_json.c_str(), mate_two_json.size());
                pair<Alignment, Alignment> ret = ff.orientation_filter(uno, dos);

                // It has to say it didn't
                // trim it
                REQUIRE(ret.first.name() == "mate_one");
                REQUIRE(ret.second.name() == "mate_two");
            }

            SECTION("Discorant read filter returns empty alignments for discordant reads in inverse mode."){

            }
        }
    }
}
