/// \file allele_splicer.cpp
///
/// Unit tests for the sequences AlleleSplicer resolves for the strands of a site and for its
/// reference allele.

#include "../../src/unittest/catch.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../allele_splicer.hpp"
#include "support/decomposition_pair.hpp"
#include "support/nested_graph.hpp"
#include "support/staged_nested_sites.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

TEST_CASE("AlleleSplicer resolves each strand through the sites nested under a site",
          "[allele_splicer]") {
    StagedNestedSites fixture;
    SECTION("A homozygous reference site takes each strand's alleles at the sites under it") {
        // The top site carries its reference walk on both strands. The middle site carries the
        // reference on strand 0 and node 6 on strand 1. The inner site is on strand 0 alone,
        // where it leaves out node 4.
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.phase(fixture.inner_key, 1, -1, 1, 0);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;

        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCATGGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATTGCTGA");
        REQUIRE(splicer.strand_traversal(fixture.site(fixture.top_key), 1) == 0);

        // The nested sites resolve the same way.
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.middle_key)) == "TGCTGAGCAG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.middle_key), 0) == "TGGCAG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.middle_key), 1) == "TTG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.inner_key), 0) == "GGCA");
        // Strand 1 does not cross the inner site.
        REQUIRE(!splicer.strand_sequence(fixture.site(fixture.inner_key), 1));
        REQUIRE(!splicer.strand_traversal(fixture.site(fixture.inner_key), 1));

        // The record of the top site spells no allele differently, since both strands carry the
        // reference.
        AlleleSplicer::SplicedAlleles spliced;
        REQUIRE(!splicer.splice(fixture.site(fixture.top_key), spliced));
    }
    SECTION("A strand that skips a nested site takes nothing from it") {
        // Strand 0 of the top site goes from node 1 to node 8. The middle site is on strand 1.
        fixture.phase(fixture.top_key, 1, 0);
        fixture.phase(fixture.middle_key, 1, -1, 1, 1);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCACTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATTGCTGA");
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
        REQUIRE(!splicer.strand_sequence(fixture.site(fixture.middle_key), 0));
    }
    SECTION("A change inside a node reaches the strand that carries it") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 0);
        fixture.phase(fixture.inner_key, 0, 0);
        fixture.phase(fixture.edit_key, 1, 0);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCATGCAGAGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATGCTGAGCAGCTGA");
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
    }
    SECTION("A site with no phase has no strands") {
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.index();
        REQUIRE(!fixture.splicer.strand_sequence(fixture.site(fixture.top_key), 0));
        REQUIRE(fixture.splicer.reference_sequence(fixture.site(fixture.top_key))
                == "GCATGCTGAGCAGCTGA");
    }
}

}
}
