#include <iostream>
#include "../algorithms/rmq.hpp"
#include "catch.hpp"
#include <vector>
#include <algorithm>

namespace vg {
namespace unittest {

using namespace vg::algorithms;

TEST_CASE("ReadCoordinateRMQAVL basic operations", "[rmq]") {
    // Basic RMQ with gap_extend = 0 (standard RMQ behavior)
    ReadCoordinateRMQAVL rmq(0);

    SECTION("Empty RMQ returns nothing") {
        auto results = rmq.query_top_k_candidates(100, 50, 5);
        REQUIRE(results.empty());
    }

    SECTION("Single item query") {
        rmq.insert(10, 100, 1);

        auto results = rmq.query_top_k_candidates(15, 10, 5);
        REQUIRE(results.size() == 1);
        REQUIRE(results[0] == 1);

        auto out_of_range = rmq.query_top_k_candidates(25, 10, 5);
        REQUIRE(out_of_range.empty());
    }

    SECTION("Multiple items top-k") {
        rmq.insert(10, 50, 1);
        rmq.insert(20, 100, 2);
        rmq.insert(30, 80, 3);
        rmq.insert(40, 120, 4);
        rmq.insert(50, 90, 5);

        // All items in range [0, 100], asking for top 3
        auto results = rmq.query_top_k_candidates(100, 100, 3);
        REQUIRE(results.size() == 3);
        // Scores: 4(120), 2(100), 5(90)
        REQUIRE(results[0] == 4);
        REQUIRE(results[1] == 2);
        REQUIRE(results[2] == 5);
    }

    SECTION("Range constraints") {
        rmq.insert(10, 100, 1);
        rmq.insert(20, 200, 2);
        rmq.insert(30, 300, 3);
        rmq.insert(40, 400, 4);

        // Range [15, 35], asking for top 2
        auto results = rmq.query_top_k_candidates(35, 20, 2);
        REQUIRE(results.size() == 2);
        // Range items: 20(200), 30(300). Top 2: 30, 20
        REQUIRE(results[0] == 3);
        REQUIRE(results[1] == 2);
    }
}

TEST_CASE("ReadCoordinateRMQAVL A&O logic (Separable Penalty)", "[rmq]") {
    // Gap cost 1 per bp
    ReadCoordinateRMQAVL rmq(1);

    SECTION("Closer item beats better but further item") {
        // Candidate A: coord 10, score 100.
        // Candidate B: coord 50, score 95.
        // Querying at coord 60.
        // Penalized Score A = 100 - (60-10)*1 = 100 - 50 = 50
        // Penalized Score B = 95 - (60-50)*1 = 95 - 10 = 85
        // So B should win even if its raw score is 95 < 100.

        // formula = score + gap_extend * coord
        // Val A = 100 + 10 * 1 = 110
        // Val B = 95 + 50 * 1 = 145
        // Query at 60 maximizes (Val - 60 * 1).

        rmq.insert(10, 100, 1);
        rmq.insert(50, 95, 2);

        auto results = rmq.query_top_k_candidates(60, 60, 1);
        REQUIRE(results.size() == 1);
        REQUIRE(results[0] == 2); // B should win
    }
}

TEST_CASE("ReadCoordinateRMQAVL balance and pruning", "[rmq]") {
    ReadCoordinateRMQAVL rmq(0);
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        rmq.insert(i, i, i);
    }

    auto results = rmq.query_top_k_candidates(1000, 1000, 10);
    REQUIRE(results.size() == 10);
    for (int i = 0; i < 10; ++i) {
        REQUIRE(results[i] == (size_t)(999 - i));
    }
}

}
}
