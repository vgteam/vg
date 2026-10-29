/// \file child_crossing_finder.cpp
///
/// Unit tests for ChildCrossingFinder, on the snarls and edits of an overlay.

#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../child_crossing_finder.hpp"
#include "support/decomposition_pair.hpp"
#include "../edit_snarl_decomposition_overlay.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "support/nested_graph.hpp"

namespace vg {
namespace unittest {

using namespace std;
using multipass::ChildCrossing;
using multipass::ChildCrossingFinder;

namespace {

/// The nested graph, its decomposition, and an overlay of edits laid over it.
struct NestedSites {
    NestedGraph nested;
    IntegratedSnarlFinder finder;
    DecompositionPair pair;
    EditSnarlDecompositionOverlay overlay;

    NestedSites() : finder(nested.graph), pair(nested.graph, finder),
                    overlay(pair.adapter, nested.graph) {}

    /// Make a walk of the given node IDs, each read forward, or backward where negative.
    Traversal walk(const vector<int>& ids) const {
        Traversal out;
        for (int id : ids) {
            out.push_back(nested.handle(abs(id), id < 0));
        }
        return out;
    }

    /// Get the entry steps of `crossings`.
    static vector<int> steps(const vector<ChildCrossing>& crossings) {
        vector<int> out;
        for (const ChildCrossing& c : crossings) {
            out.push_back(c.entry_step);
        }
        return out;
    }
};

}

TEST_CASE("A walk crosses a child snarl by visiting one bound and then the other",
          "[child_crossing_finder]") {
    NestedSites sites;
    const net_handle_t middle = find_snarl(sites.overlay, sites.nested.graph, 2, 7);
    const net_handle_t inner = find_snarl(sites.overlay, sites.nested.graph, 3, 5);

    SECTION("A walk through the snarl crosses it once, entering at its first bound") {
        const Traversal walk = sites.walk({1, 2, 3, 5, 7, 8});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        const vector<ChildCrossing> crossed = finder.crossings(middle);
        REQUIRE(NestedSites::steps(crossed) == vector<int>{1});
        // Node 1 is 3 bases long.
        REQUIRE(crossed.front().entry_bases == 3);
        // Either orientation of the snarl is the same child.
        REQUIRE(NestedSites::steps(finder.crossings(sites.overlay.flip(middle))) == vector<int>{1});
        REQUIRE(NestedSites::steps(finder.crossings(inner)) == vector<int>{2});
    }
    SECTION("A walk that skips the snarl does not cross it") {
        const Traversal walk = sites.walk({1, 8});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        REQUIRE(finder.crossings(middle).empty());
    }
    SECTION("A walk read backward crosses the snarl too") {
        const Traversal walk = sites.walk({-8, -7, -6, -2, -1});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        REQUIRE(NestedSites::steps(finder.crossings(middle)) == vector<int>{1});
        REQUIRE(finder.crossings(inner).empty());
    }
    SECTION("A walk that passes through the snarl twice crosses it twice") {
        const Traversal walk = sites.walk({2, 3, 5, 7, 2, 6, 7});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        const vector<ChildCrossing> crossed = finder.crossings(middle);
        REQUIRE(NestedSites::steps(crossed) == vector<int>{0, 4});
        // T, G, GCA and G come before the second entry.
        REQUIRE(crossed[1].entry_bases == 6);
    }
}

TEST_CASE("A walk crosses an edit inside a node at each visit to the node",
          "[child_crossing_finder]") {
    NestedSites sites;
    // Node 4 is CTGA. The edit starts after its first base.
    const net_handle_t edit =
        sites.overlay.add_edit(sites.nested.handle(4), 1, sites.nested.handle(4), 2);

    SECTION("A visit read forward enters the edit at its first base") {
        const Traversal walk = sites.walk({3, 4, 5});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        const vector<ChildCrossing> crossed = finder.crossings(edit);
        REQUIRE(NestedSites::steps(crossed) == vector<int>{1});
        // G, then C.
        REQUIRE(crossed.front().entry_bases == 2);
    }
    SECTION("A visit read backward counts the node's bases from the other end") {
        const Traversal walk = sites.walk({-5, -4, -3});
        const ChildCrossingFinder finder(sites.overlay, sites.nested.graph, walk);
        const vector<ChildCrossing> crossed = finder.crossings(edit);
        REQUIRE(NestedSites::steps(crossed) == vector<int>{1});
        // GCA, then the two bases of node 4 after the edit's first base.
        REQUIRE(crossed.front().entry_bases == 5);
    }
    SECTION("Two visits cross it twice, and a walk that misses the node not at all") {
        const Traversal twice = sites.walk({3, 4, 5, 4});
        REQUIRE(NestedSites::steps(
                    ChildCrossingFinder(sites.overlay, sites.nested.graph, twice).crossings(edit))
                == vector<int>{1, 3});
        const Traversal missing = sites.walk({3, 5});
        const ChildCrossingFinder missing_finder(sites.overlay, sites.nested.graph, missing);
        REQUIRE(missing_finder.crossings(edit).empty());
    }
}

TEST_CASE("The crossing mask and the entry along a genotype read the walks that cross a child",
          "[child_crossing_finder]") {
    NestedSites sites;
    const net_handle_t edit =
        sites.overlay.add_edit(sites.nested.handle(4), 1, sites.nested.handle(4), 2);
    const net_handle_t inner = find_snarl(sites.overlay, sites.nested.graph, 3, 5);
    const vector<Traversal> travs{sites.walk({3, 5}), sites.walk({3, 4, 5}),
                                  sites.walk({3, 4, 5, 4})};
    vector<ChildCrossingFinder> walks;
    for (const Traversal& t : travs) {
        walks.emplace_back(sites.overlay, sites.nested.graph, t);
    }

    SECTION("The mask sets the bit of each walk that crosses the child") {
        bool known = false;
        REQUIRE(multipass::crossing_mask(walks, edit, &known) == 0b110);
        REQUIRE(known);
        REQUIRE(multipass::crossing_mask(walks, inner, &known) == 0b111);
    }
    SECTION("More than 64 walks make the mask unknown") {
        const vector<ChildCrossingFinder> many(65, walks.front());
        bool known = true;
        REQUIRE(multipass::crossing_mask(many, inner, &known) == 0);
        REQUIRE_FALSE(known);
    }
    SECTION("The entry is along the first walk of the genotype that crosses the child") {
        unordered_map<const Traversal*, ChildCrossingFinder> finders;
        REQUIRE(multipass::entry_bases_along_genotype(sites.overlay, sites.nested.graph, travs,
                                                      {0, 1}, edit, finders)
                == 2);
        REQUIRE(multipass::entry_bases_along_genotype(sites.overlay, sites.nested.graph, travs,
                                                      {0}, edit, finders)
                == 0);
    }
}

}
}
