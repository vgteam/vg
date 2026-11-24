/// \file staged_site.cpp
///
/// Unit tests for StagedSiteTable.

#include "../../src/unittest/catch.hpp"
#include "../staged_site.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

TEST_CASE("StagedSiteTable::hand_off frees the panel alleles of the sites it holds back",
          "[staged_site]") {
    StagedSiteTable table;
    table.start(1);

    StagedSite site;
    site.placement.reported_inline = true;
    site.panel_cache = vector<int>{0, 1, -1, 1};
    table.add_nested(std::move(site));
    table.gather_nested();

    StagedSiteTable::HandOff counts = table.hand_off();
    REQUIRE(counts.inline_unrendered == 1);
    REQUIRE(table.unrendered().size() == 1);

    const StagedSite& held = table.unrendered().front();
    REQUIRE(!held.panel_cache);
}

}
}
