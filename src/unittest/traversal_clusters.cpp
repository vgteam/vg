/// \file traversal_clusters.cpp
///
/// Unit tests for traversal similarity, the metric behind -L/--cluster in vg deconstruct,
/// vg call and vg simplify.
///
/// These exercise weighted_traversal_similarity directly rather than through a caller.  That
/// matters: the pure-deletion blindness this file's second half pins went unnoticed through
/// several reviews precisely because the metric was only ever reached through two layers of
/// caller, where a similarity of 0 is indistinguishable from "these alleles really are unrelated".
///

#include <iostream>
#include <set>
#include "../traversal_clusters.hpp"
#include "catch.hpp"
#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("traversal similarity handles pure deletions", "[jaccard][traversal_clusters]") {

    // A 60bp site: node 2 is 59bp, node 3 is 1bp, with 10bp boundaries.
    //   ref    = 1 -> 2 -> 3 -> 4     interior {2,3}, 60bp
    //   del59  = 1 -> 3 -> 4          interior {3},    1bp  (deletes the 59bp node)
    //   del60  = 1 -> 4               interior {},     0bp  (deletes the whole site)
    //   sub59  = 1 -> 5 -> 3 -> 4     interior {5,3}, 60bp  (59bp substitution)
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle(string(10, 'A'));
    handle_t h2 = graph.create_handle(string(59, 'C'));
    handle_t h3 = graph.create_handle("G");
    handle_t h4 = graph.create_handle(string(10, 'T'));
    handle_t h5 = graph.create_handle(string(59, 'A'));

    multiset<handle_t> ref{h2, h3};
    multiset<handle_t> del59{h3};
    multiset<handle_t> del60{};
    multiset<handle_t> sub59{h5, h3};
    multiset<handle_t> no_site{};

    auto sim = [&](const multiset<handle_t>& a, const multiset<handle_t>& b,
                   const multiset<handle_t>& site) {
        return weighted_traversal_similarity(&graph, a, b, site);
    };

    SECTION("identical alleles score exactly 1, at any site") {
        REQUIRE(sim(ref, ref, ref) == Approx(1.0));
        REQUIRE(sim(del59, del59, ref) == Approx(1.0));
        // the 0/0 input: two pure deletions of the same site.  1, not NaN -- a NaN fails both the
        // > and the >= comparisons in cluster_traversals and would make the traversal permanently
        // unclusterable, which is the same class of bug this metric change exists to fix.
        REQUIRE(sim(del60, del60, ref) == Approx(1.0));
        REQUIRE(sim(del60, del60, no_site) == Approx(1.0));
    }

    SECTION("two alleles that both carry sequence ignore the site entirely") {
        // regression guard: this is the pre-existing behaviour and must not drift, or every
        // existing -L threshold in the wild silently changes meaning
        for (const multiset<handle_t>* site : {&no_site, &ref, &sub59}) {
            REQUIRE(sim(ref, sub59, *site) == Approx(1.0 / 119.0));
            REQUIRE(sim(ref, del59, *site) == Approx(1.0 / 60.0));
        }
    }

    SECTION("THE DEFECT: a pure deletion is similar to a near-identical deletion") {
        // 59 of the site's 60bp are deleted by both.  Before the fix this was 0 at every
        // threshold, because a 2-visit traversal kept its boundary handles and was therefore
        // disjoint from every other allele by construction.
        REQUIRE(sim(del59, del60, ref) == Approx(59.0 / 60.0));
        REQUIRE(sim(del60, del59, ref) == Approx(59.0 / 60.0));
    }

    SECTION("a pure deletion stays apart from alleles that fill the site") {
        // Both of these carry a full site's worth of sequence while the deletion carries none, so
        // the fraction of the site they agree on deleting is 0.  Only the amount DELETED is
        // compared -- a 59bp substitution is not "59bp similar" to a 59bp deletion just because
        // the two numbers match.
        REQUIRE(sim(ref, del60, ref) == Approx(0.0));
        REQUIRE(sim(sub59, del60, ref) == Approx(0.0));
    }

    SECTION("without a site a pure deletion falls back to pairwise scoring") {
        // i.e. the pre-fix answer, which is what callers that cannot supply a reference get
        REQUIRE(sim(del59, del60, no_site) == Approx(0.0));
    }

    SECTION("the metric is symmetric and bounded") {
        vector<const multiset<handle_t>*> all{&ref, &del59, &del60, &sub59};
        for (const multiset<handle_t>* a : all) {
            for (const multiset<handle_t>* b : all) {
                double ab = sim(*a, *b, ref);
                double ba = sim(*b, *a, ref);
                REQUIRE(ab == Approx(ba));
                REQUIRE(ab >= 0.0);
                REQUIRE(ab <= 1.0);
            }
        }
    }

    SECTION("orientation matters: an inversion shares no handles with its forward allele") {
        multiset<handle_t> fwd{h2};
        multiset<handle_t> rev{graph.flip(h2)};
        REQUIRE(sim(fwd, rev, ref) == Approx(0.0));
        // and the site does not rescue it -- neither allele is a pure deletion
        multiset<handle_t> big_site{h2, h3, h5};
        REQUIRE(sim(fwd, rev, big_site) == Approx(0.0));
    }

    SECTION("multiset semantics distinguish cycle copy number") {
        multiset<handle_t> once{h3};
        multiset<handle_t> twice{h3, h3};
        REQUIRE(sim(once, twice, ref) == Approx(0.5));
    }
}

}
}
