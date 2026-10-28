/// \file round_history.cpp
/// Tests for RoundHistory, the record of the states that re-genotyping rounds reach.

#include "../../src/unittest/catch.hpp"
#include "../round_history.hpp"

namespace vg {
namespace unittest {

using namespace std;
using multipass::RoundHistory;

TEST_CASE("RoundHistory finds the first round that reached a state", "[round_history]") {
    RoundHistory history;
    RoundHistory::State one{{10, {0, 1, 2}}, {20, {1, 1, 2}}};
    RoundHistory::State two{{10, {0, 0, 2}}, {20, {1, 1, 2}}};
    history.remember(one);
    history.remember(two);
    REQUIRE(history.rounds() == 2);
    REQUIRE(history.first_round_with(one) == 1);
    REQUIRE(history.first_round_with(two) == 2);
    REQUIRE(history.first_round_with(RoundHistory::State{{10, {0, 1, 2}}}) == 0);

    // The same state filed in another order is the same state.
    RoundHistory::State reordered;
    reordered[20] = {1, 1, 2};
    reordered[10] = {0, 0, 2};
    REQUIRE(history.first_round_with(reordered) == 2);
}

TEST_CASE("RoundHistory does not match different states with equal hashes",
          "[round_history]") {
    // Mixing a site's fields as base-1000003 digits makes these two sites hash alike, since
    // 1 * 1000003 + 0 == 0 * 1000003 + 1000003.
    RoundHistory history;
    history.remember(RoundHistory::State{{0, {-3, -2, 0}}});
    REQUIRE(history.first_round_with(RoundHistory::State{{0, {-3, -3, 1000003}}}) == 0);
}

}
}
