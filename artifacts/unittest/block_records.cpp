/// \file block_records.cpp
/// Unit tests for writing a site as block records.

#include <stdexcept>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>

#include "../../src/unittest/catch.hpp"
#include "../block_records.hpp"

namespace vg {
namespace unittest {

using namespace std;
using namespace vg::multipass;

TEST_CASE("BlockRecordWriter rejects a genotype entry that names no traversal",
          "[block_records]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A", 1);
    handle_t n2 = graph.create_handle("C", 2);
    handle_t n3 = graph.create_handle("G", 3);
    handle_t n4 = graph.create_handle("T", 4);
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    bdsg::PositionOverlay positioned(&graph);

    vector<Traversal> travs = {{n1, n2, n4}, {n1, n3, n4}};
    SiteChildren children;
    children.known = true;
    SiteRecord record{};
    record.alleles_merged = false;

    BlockRecordWriter writer;
    writer.set_nested(true);
    writer.set_enabled(true);
    // Entry 5 is past the two traversals, while entry 1 differs from the reference, so the site
    // has a difference block to cut.
    REQUIRE_THROWS_AS(writer.write(positioned, children, travs, {1, 5}, 0, "s", nullptr, record,
                                   GLLayout::IMajor, false,
                                   [](vcflib::Variant&, size_t) { return true; }),
                      std::out_of_range);
}

}
}
