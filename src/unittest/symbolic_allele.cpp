/// \file symbolic_allele.cpp
/// Unit tests for projecting a traversal into symbolic form.
///
/// The cases that matter are the ones where getting it wrong is silent. Collapsing too eagerly makes
/// two genuinely different alleles compare equal, so a real variant is dropped with no record and no
/// warning. Collapsing too little leaves the nested variation baked into a long allele, which is the
/// behaviour this exists to remove. Both directions are tested.

#include <vector>

#include "catch.hpp"
#include "../symbolic_allele.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A snarl from start to end, forward on both boundaries.
static Snarl make_snarl(nid_t start, nid_t end) {
    Snarl s;
    s.mutable_start()->set_node_id(start);
    s.mutable_end()->set_node_id(end);
    return s;
}

/// A traversal over the given node ids, all forward.
static SnarlTraversal make_trav(const vector<nid_t>& nodes) {
    SnarlTraversal t;
    for (nid_t n : nodes) {
        t.add_visit()->set_node_id(n);
    }
    return t;
}

/// A SnarlManager holding one top-level snarl and the given children of it.
static unique_ptr<SnarlManager> make_manager(const Snarl& top, const vector<Snarl>& children) {
    vector<Snarl> all;
    all.push_back(top);
    for (Snarl c : children) {
        // The manager builds the hierarchy from the parent links, so they have to be set.
        *c.mutable_parent()->mutable_start() = top.start();
        *c.mutable_parent()->mutable_end() = top.end();
        all.push_back(c);
    }
    return unique_ptr<SnarlManager>(new SnarlManager(all.begin(), all.end()));
}

TEST_CASE("A traversal with no child chain is its own symbolic form", "[symbolic_allele]") {
    Snarl top = make_snarl(1, 5);
    auto mgr = make_manager(top, {});
    SnarlTraversal t = make_trav({1, 2, 3, 5});

    SymbolicAllele a = symbolic_allele(t, top, *mgr);
    REQUIRE(a.size() == 4);
    for (const SymbolicStep& s : a) {
        REQUIRE_FALSE(s.is_chain());
    }
    REQUIRE_FALSE(has_child_chain(t, top, *mgr));
}

TEST_CASE("Traversals differing only inside a child chain are symbolically equal",
          "[symbolic_allele]") {
    // This is the whole point. 1 -> [child 2..4] -> 5, and the two traversals take different routes
    // through the child. On the real data these are the 4,710 bp alleles differing at three bases.
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    SnarlTraversal via_a = make_trav({1, 2, 3, 4, 5});
    SnarlTraversal via_b = make_trav({1, 2, 30, 4, 5});

    REQUIRE(has_child_chain(via_a, top, *mgr));
    REQUIRE(symbolically_equal(via_a, via_b, top, *mgr));

    // The form is [1, chain(2..4), 4, 5]: the exit boundary is emitted as a plain node because it
    // belongs to whatever follows the chain as much as to the chain, so it is not consumed. Both
    // traversals produce that same form, which is what makes them equal.
    SymbolicAllele a = symbolic_allele(via_a, top, *mgr);
    REQUIRE(a.size() == 4);
    REQUIRE(a[1].is_chain());
    REQUIRE(a[1].id == 2);
    REQUIRE(a[1].end_id == 4);
    REQUIRE_FALSE(a[2].is_chain());
    REQUIRE(a[2].id == 4);
}

TEST_CASE("A traversal that skips a child chain stays symbolically distinct", "[symbolic_allele]") {
    // A real deletion must survive: it does not cross the chain at all, so it must not collapse
    // into the reference. If this ever fails, deletions vanish silently.
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    SnarlTraversal crosses = make_trav({1, 2, 3, 4, 5});
    SnarlTraversal skips = make_trav({1, 5});

    REQUIRE_FALSE(symbolically_equal(crosses, skips, top, *mgr));
    REQUIRE_FALSE(has_child_chain(skips, top, *mgr));
}

TEST_CASE("The site's own boundaries are never collapsed", "[symbolic_allele]") {
    // A traversal enters at the site's own start. Symbolising that would reduce every allele to a
    // single symbol and make all of them equal -- every variant in the graph would disappear.
    Snarl top = make_snarl(1, 5);
    auto mgr = make_manager(top, {});

    SnarlTraversal a = make_trav({1, 2, 5});
    SnarlTraversal b = make_trav({1, 3, 5});

    REQUIRE_FALSE(symbolically_equal(a, b, top, *mgr));
}

TEST_CASE("A chain entered and never left is not collapsed", "[symbolic_allele]") {
    // Truncated or cyclic traversals happen. Swallowing the tail would drop everything after the
    // entry point and make unrelated alleles compare equal, so such a step stays a plain node.
    Snarl top = make_snarl(1, 9);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    SnarlTraversal truncated = make_trav({1, 2, 3, 9});   // enters the child, never reaches node 4
    SymbolicAllele a = symbolic_allele(truncated, top, *mgr);
    for (const SymbolicStep& s : a) {
        REQUIRE_FALSE(s.is_chain());
    }
    REQUIRE(a.size() == 4);
}

TEST_CASE("Two different child chains give different symbols", "[symbolic_allele]") {
    Snarl top = make_snarl(1, 9);
    Snarl first = make_snarl(2, 4);
    Snarl second = make_snarl(5, 7);
    auto mgr = make_manager(top, {first, second});

    SnarlTraversal both = make_trav({1, 2, 3, 4, 5, 6, 7, 9});
    SymbolicAllele a = symbolic_allele(both, top, *mgr);

    vector<SymbolicStep> chains;
    for (const SymbolicStep& s : a) {
        if (s.is_chain()) {
            chains.push_back(s);
        }
    }
    REQUIRE(chains.size() == 2);
    REQUIRE(chains[0] != chains[1]);
}

TEST_CASE("A visit already carrying a snarl is taken as a symbol", "[symbolic_allele]") {
    // The protobuf allows a Visit to hold a Snarl instead of a node, and the deprecated
    // NestedFlowCaller emitted them, so the encoder has to accept traversals already in that form.
    Snarl top = make_snarl(1, 9);
    auto mgr = make_manager(top, {});

    SnarlTraversal t;
    t.add_visit()->set_node_id(1);
    Visit* v = t.add_visit();
    v->mutable_snarl()->mutable_start()->set_node_id(2);
    v->mutable_snarl()->mutable_end()->set_node_id(4);
    t.add_visit()->set_node_id(9);

    SymbolicAllele a = symbolic_allele(t, top, *mgr);
    REQUIRE(a.size() == 3);
    REQUIRE(a[1].is_chain());
    REQUIRE(a[1].id == 2);
    REQUIRE(a[1].end_id == 4);
}

TEST_CASE("Symbolic alleles hash by value", "[symbolic_allele]") {
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    SnarlTraversal via_a = make_trav({1, 2, 3, 4, 5});
    SnarlTraversal via_b = make_trav({1, 2, 30, 4, 5});

    hash<SymbolicAllele> h;
    REQUIRE(h(symbolic_allele(via_a, top, *mgr)) == h(symbolic_allele(via_b, top, *mgr)));
}

}
}
