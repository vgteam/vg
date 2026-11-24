/// \file weakly_connected_components_dense.cpp
///
/// Unit tests for weakly_connected_components_dense.

#include <unordered_set>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <handlegraph/algorithms/weakly_connected_components.hpp>

#include "../algorithms/weakly_connected_components_dense.hpp"
#include "catch.hpp"
#include "support/random_graph.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("weakly_connected_components_dense finds the components handlealgs does, in its order",
          "[weakly_connected_components_dense]") {
    for (int64_t seed = 0; seed < 20; ++seed) {
        bdsg::HashGraph graph;
        // Several components, of different sizes.
        random_graph({100 + seed, 40, 250, 7}, 5, 8, &graph);
        INFO("seed " << seed);
        REQUIRE(algorithms::weakly_connected_components_dense(graph)
                == handlealgs::weakly_connected_components(&graph));
    }
}

TEST_CASE("weakly_connected_components_dense gives an empty graph no components",
          "[weakly_connected_components_dense]") {
    bdsg::HashGraph graph;
    REQUIRE(algorithms::weakly_connected_components_dense(graph).empty());
}

}
}
