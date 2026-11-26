/// \file staged_site.cpp
///
/// Unit tests for StagedSiteTable.

#include <algorithm>

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

TEST_CASE("StagedSiteTable::remove_with_nested removes the named sites and the sites under them",
          "[staged_site]") {
    StagedSiteTable table;
    table.start(1);
    auto site = [](size_t key, size_t parent) {
        StagedSite made;
        made.record_key = key;
        made.placement.parent_record_key = parent;
        made.placement.level = parent == 0 ? 0 : 1;
        return made;
    };
    // Two top-level sites, 1 and 2. Under 1 are 11 and 12, under 11 is 111, and under 2 is 21.
    table.add_top_level(site(1, 0));
    table.add_top_level(site(2, 0));
    table.add_nested(site(11, 1));
    table.add_nested(site(21, 2));
    table.add_nested(site(12, 1));
    table.add_nested(site(111, 11));
    table.gather_nested();

    SECTION("A top-level site goes with everything under it") {
        vector<size_t> gone;
        const size_t removed = table.remove_with_nested(
            {1}, [&](const StagedSite& site) { gone.push_back(site.record_key); });
        REQUIRE(removed == 4);
        std::ranges::sort(gone);
        REQUIRE(gone == vector<size_t>{1, 11, 12, 111});
        REQUIRE(table.queued_count() == 1);
        REQUIRE(table.nested().size() == 1);
        REQUIRE(table.children_of(1) == nullptr);
        // The sites that stay are indexed where they are now.
        REQUIRE(table.children_of(2) != nullptr);
        REQUIRE(*table.children_of(2) == vector<size_t>{0});
        REQUIRE(table.nested()[0].record_key == 21);
    }

    SECTION("A key that is not staged still names the sites under it") {
        // As where a site made of a run has taken the place of site 11, under another key.
        table.nested()[0].record_key = 13;
        const size_t removed = table.remove_with_nested({11}, [](const StagedSite&) {
            // Nothing to note.
        });
        REQUIRE(removed == 1);
        REQUIRE(table.nested().size() == 3);
        REQUIRE(table.children_of(11) == nullptr);
        REQUIRE(table.children_of(1)->size() == 2);
    }
}

}
}
