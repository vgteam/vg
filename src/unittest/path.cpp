/// \file path.cpp
///  
/// unit tests for Paths and their utility functions
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/reference_path_overlay.hpp>
#include "../path.hpp"
#include "../vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Path simplification tolerates adjacent insertions and deletions", "[path]") {

    string path_string = R"(
        {
            "mapping": [
                {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "68"}},
                {"edit": [{"sequence": "AAGG", "to_length": 4}, {"from_length": 3}], "position": {"node_id": "67"}},
                {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "66"}}
            ]
        }
    )";

    Path path;
    json2pb(path, path_string.c_str(), path_string.size());

    // Simplify without replacing deletions with skips
    auto simple = simplify(path, false);
    
    // We need to still touch all the nodes after simplification.
    REQUIRE(simple.mapping_size() == 3);
    REQUIRE(simple.mapping(0).position().node_id() == 68);
    REQUIRE(simple.mapping(1).position().node_id() == 67);
    REQUIRE(simple.mapping(2).position().node_id() == 66);

}

TEST_CASE("for_each_overlapping_subpath() computes overlapping ranges correctly", "[path]") {

    bdsg::HashGraph base_graph;
    handle_t base_h1 = base_graph.create_handle("GAT");
    handle_t base_h2 = base_graph.create_handle("T");
    handle_t base_h3 = base_graph.create_handle("ACA");

    SECTION("With a full base path") {
        path_handle_t base_ph1 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE);
        base_graph.append_step(base_ph1, base_h1);
        base_graph.append_step(base_ph1, base_h2);
        base_graph.append_step(base_ph1, base_h3);

        bdsg::ReferencePathOverlay ref_graph(&base_graph, {});

        REQUIRE(ref_graph.has_path("GRCh38#0#chr1"));
        path_handle_t ref_ph1 = ref_graph.get_path_handle("GRCh38#0#chr1");
        
        SECTION("A subregion finds the base path") {

            // When we look up GRCh38#0#chr1:1-3, 0-based, end-inclusive, we shoud see that path and the right end-exclusive range on it.
            Region target_region {"GRCh38#0#chr1", 1, 3};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                REQUIRE(path == ref_ph1);
                REQUIRE(start_offset == 1);
                REQUIRE(past_end_offset == 4);
                results++;
                return true;
            });
            REQUIRE(results == 1);
        }

        SECTION("An open left end finds the base path") {

            Region target_region {"GRCh38#0#chr1", -1, 3};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                REQUIRE(path == ref_ph1);
                REQUIRE(start_offset == 0);
                REQUIRE(past_end_offset == 4);
                results++;
                return true;
            });
            REQUIRE(results == 1);
        }

        SECTION("An open right end finds the base path") {

            Region target_region {"GRCh38#0#chr1", 1, -1};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                REQUIRE(path == ref_ph1);
                REQUIRE(start_offset == 1);
                REQUIRE(past_end_offset == 7);
                results++;
                return true;
            });
            REQUIRE(results == 1);
        }

        SECTION("No bounds finds the base path") {

            Region target_region {"GRCh38#0#chr1", -1, -1};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                REQUIRE(path == ref_ph1);
                REQUIRE(start_offset == 0);
                REQUIRE(past_end_offset == 7);
                results++;
                return true;
            });
            REQUIRE(results == 1);
        }
    }

    SECTION("With parts of a base path") {
        path_handle_t base_ph1 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {0, 3});
        base_graph.append_step(base_ph1, base_h1);
        path_handle_t base_ph2 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {3, 4});
        base_graph.append_step(base_ph2, base_h2);
        path_handle_t base_ph3 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {4, 7});
        base_graph.append_step(base_ph3, base_h3);

        bdsg::ReferencePathOverlay ref_graph(&base_graph, {});
        
        // TODO: GBZ will name this differently, as just "GRCh38#0#chr1" and without the subrange on it.
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[0-3]"));
        path_handle_t ref_ph1 = ref_graph.get_path_handle("GRCh38#0#chr1[0-3]");
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[3-4]"));
        path_handle_t ref_ph2 = ref_graph.get_path_handle("GRCh38#0#chr1[3-4]");
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[4-7]"));
        path_handle_t ref_ph3 = ref_graph.get_path_handle("GRCh38#0#chr1[4-7]");
        
        SECTION("A subregion of the base path finds the right subpaths") {

            // When we look up GRCh38#0#chr1:1-3, 0-based, end-inclusive, we shoud see that path and the right end-exclusive range on it.
            Region target_region {"GRCh38#0#chr1", 1, 3};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                if (path == ref_ph1) {
                    REQUIRE(start_offset == 1);
                    REQUIRE(past_end_offset == 3);
                } else if (path == ref_ph2) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 1);
                } else {
                    REQUIRE(false);
                }
                results++;
                return true;
            });
            REQUIRE(results == 2);
        }

        SECTION("An open left end finds the right subpaths") {

            Region target_region {"GRCh38#0#chr1", -1, 3};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                if (path == ref_ph1) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 3);
                } else if (path == ref_ph2) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 1);
                } else {
                    REQUIRE(false);
                }
                results++;
                return true;
            });
            REQUIRE(results == 2);
        }

        SECTION("An open right end finds the right subpaths") {

            Region target_region {"GRCh38#0#chr1", 1, -1};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                if (path == ref_ph1) {
                    REQUIRE(start_offset == 1);
                    REQUIRE(past_end_offset == 3);
                } else if (path == ref_ph2) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 1);
                } else if (path == ref_ph3) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 3);
                } else {
                    REQUIRE(false);
                }
                results++;
                return true;
            });
            REQUIRE(results == 3);
        }

        SECTION("No bounds finds all the subpaths") {

            Region target_region {"GRCh38#0#chr1", -1, -1};
            size_t results = 0;
            vg::for_each_overlapping_subpath(ref_graph, target_region, [&](const path_handle_t& path, size_t start_offset, size_t past_end_offset) {
                if (path == ref_ph1) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 3);
                } else if (path == ref_ph2) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 1);
                } else if (path == ref_ph3) {
                    REQUIRE(start_offset == 0);
                    REQUIRE(past_end_offset == 3);
                } else {
                    REQUIRE(false);
                }
                results++;
                return true;
            });
            REQUIRE(results == 3);
        }
    }
}

TEST_CASE("find_containing_subpath() works", "[path]") {
    bdsg::HashGraph base_graph;
    handle_t base_h1 = base_graph.create_handle("GAT");
    handle_t base_h2 = base_graph.create_handle("T");
    handle_t base_h3 = base_graph.create_handle("ACA");

    SECTION("With a full base path") {
        path_handle_t base_ph1 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE);
        base_graph.append_step(base_ph1, base_h1);
        base_graph.append_step(base_ph1, base_h2);
        base_graph.append_step(base_ph1, base_h3);

        bdsg::ReferencePathOverlay ref_graph(&base_graph, {});

        REQUIRE(ref_graph.has_path("GRCh38#0#chr1"));
        path_handle_t ref_ph1 = ref_graph.get_path_handle("GRCh38#0#chr1");

        path_handle_t path;
        
        SECTION("It finds the base path and populates the region") {
            Region target_region {"GRCh38#0#chr1", -1, -1};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 0);
            REQUIRE(target_region.end == 6);
            REQUIRE(path == ref_ph1);
        }

        SECTION("It finds the base path and populates the region start") {
            Region target_region {"GRCh38#0#chr1", -1, 5};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 0);
            REQUIRE(target_region.end == 5);
            REQUIRE(path == ref_ph1);
        }

        SECTION("It finds the base path and populates the region end") {
            Region target_region {"GRCh38#0#chr1", 1, -1};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 1);
            REQUIRE(target_region.end == 6);
            REQUIRE(path == ref_ph1);
        }

        SECTION("It finds the base path when asked about a fixed region") {
            Region target_region {"GRCh38#0#chr1", 3, 4};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 3);
            REQUIRE(target_region.end == 4);
            REQUIRE(path == ref_ph1);
        }

        SECTION("It finds the base path when asked about an empty region") {
            Region target_region {"GRCh38#0#chr1", 3, 2};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 3);
            REQUIRE(target_region.end == 2);
            REQUIRE(path == ref_ph1);
        }

        SECTION("It does not find the base path when asked about a region outside it") {
            Region target_region {"GRCh38#0#chr1", 100, 103};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == false);
            REQUIRE(target_region.start == 100);
            REQUIRE(target_region.end == 103);
        }
    }

    
    SECTION("With parts of a base path") {
        path_handle_t base_ph1 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {0, 3});
        base_graph.append_step(base_ph1, base_h1);
        path_handle_t base_ph2 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {3, 4});
        base_graph.append_step(base_ph2, base_h2);
        path_handle_t base_ph3 = base_graph.create_path(PathSense::REFERENCE, "GRCh38", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, {4, 7});
        base_graph.append_step(base_ph3, base_h3);

        bdsg::ReferencePathOverlay ref_graph(&base_graph, {});
        
        // TODO: GBZ will name this differently, as just "GRCh38#0#chr1" and without the subrange on it.
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[0-3]"));
        path_handle_t ref_ph1 = ref_graph.get_path_handle("GRCh38#0#chr1[0-3]");
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[3-4]"));
        path_handle_t ref_ph2 = ref_graph.get_path_handle("GRCh38#0#chr1[3-4]");
        REQUIRE(ref_graph.has_path("GRCh38#0#chr1[4-7]"));
        path_handle_t ref_ph3 = ref_graph.get_path_handle("GRCh38#0#chr1[4-7]");

        path_handle_t path;
        
        SECTION("It finds a subpath and populates the region") {
            Region target_region {"GRCh38#0#chr1[3-4]", -1, -1};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(target_region.start == 0);
            REQUIRE(target_region.end == 0);
            REQUIRE(path == ref_ph2);
        }

        SECTION("It does not find a subpath when the region spans a break") {
            Region target_region {"GRCh38#0#chr1", 2, 3};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == false);
        }

        SECTION("It finds a subpath when the region is in a subpath") {
            Region target_region {"GRCh38#0#chr1", 4, 5};
            bool result = find_containing_subpath(ref_graph, target_region, path);

            REQUIRE(result == true);
            REQUIRE(path == ref_ph3);
        }

    }
}

}
}
