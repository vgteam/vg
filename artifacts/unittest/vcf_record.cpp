/// \file vcf_record.cpp
/// Unit tests for placing sites on a reference path and trimming the alleles of a record.

#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>

#include "../../src/unittest/catch.hpp"
#include "../vcf_record.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// Make a path named "ref" through single-base nodes, visiting `steps` in order.
///
/// Each step is a node ID, negative to visit the node backward. Node IDs run from 1 to `nodes`.
static void make_ref_path(bdsg::HashGraph& graph, nid_t nodes, const vector<nid_t>& steps) {
    for (nid_t id = 1; id <= nodes; ++id) {
        graph.create_handle("A", id);
    }
    path_handle_t path = graph.create_path_handle("ref");
    for (nid_t step : steps) {
        graph.append_step(path, graph.get_handle(std::abs(step), step < 0));
    }
}

TEST_CASE("get_ref_interval finds a site whose end is the last step of the path",
          "[vcf_record]") {
    bdsg::HashGraph graph;
    // Visiting node 2 twice sends get_ref_interval to scan the path.
    make_ref_path(graph, 3, {2, 1, 2, 3});
    bdsg::PositionOverlay positioned(&graph);

    auto interval = get_ref_interval(positioned, positioned.get_handle(2), positioned.get_handle(3),
                                     "ref");
    REQUIRE(get<0>(interval) == 0);
    REQUIRE(get<1>(interval) == 3);
    REQUIRE_FALSE(get<2>(interval));
}

TEST_CASE("get_ref_interval scans forward for the end after a backward scan fails",
          "[vcf_record]") {
    bdsg::HashGraph graph;
    // The path first visits the start node backward, and a backward scan from there finds nothing.
    // The second visit is forward, and the end follows it.
    make_ref_path(graph, 5, {-2, 4, 2, 3, 5});
    bdsg::PositionOverlay positioned(&graph);

    auto interval = get_ref_interval(positioned, positioned.get_handle(2), positioned.get_handle(3),
                                     "ref");
    // The interval starts at the second visit to the start node, which the end follows.
    REQUIRE(get<0>(interval) == 2);
    REQUIRE(get<1>(interval) == 3);
    REQUIRE(positioned.get_handle_of_step(get<4>(interval)) == positioned.get_handle(3));
}

TEST_CASE("flatten_common_allele_ends leaves a base when len_override is longer than an allele",
          "[vcf_record]") {
    vcflib::Variant variant;
    variant.position = 10;
    variant.ref = "ACGT";
    variant.alt = {"AC"};
    variant.alleles = {"ACGT", "AC"};

    flatten_common_allele_ends(variant, false, 5);

    REQUIRE(variant.position == 11);
    REQUIRE(variant.ref == "CGT");
    REQUIRE(variant.alt == vector<string>{"C"});
}

}
}
