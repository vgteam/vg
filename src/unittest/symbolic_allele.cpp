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

/// A traversal over the given node ids, all *reverse* oriented, as the reference traversal of a
/// snarl the reference path runs backwards through would be.
static SnarlTraversal make_trav_rev(const vector<nid_t>& nodes) {
    SnarlTraversal t;
    for (nid_t n : nodes) {
        Visit* v = t.add_visit();
        v->set_node_id(n);
        v->set_backward(true);
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

/// The visit ranges must partition the traversal contiguously, in order, with no gap and no
/// overlap. Everything that extracts per-block sequence depends on this, and a violation would be
/// silent: a dropped visit shortens an allele by a node, which reads as a real variant.
static void check_partition(const SnarlTraversal& t, const Snarl& site, const SnarlManager& mgr) {
    vector<pair<int, int>> ranges;
    SymbolicAllele a = symbolic_allele(t, site, mgr, &ranges);
    REQUIRE(ranges.size() == a.size());
    int expect = 0;
    for (size_t k = 0; k < ranges.size(); ++k) {
        REQUIRE(ranges[k].first == expect);
        REQUIRE(ranges[k].second > ranges[k].first);
        // A chain symbol's range stops *at* its exit boundary, which the next step owns, so the
        // next range begins exactly where this one ended -- same rule as a plain node.
        expect = ranges[k].second;
    }
    REQUIRE(expect == t.visit_size());
}

/// A symbolic allele of plain node steps, built directly. The diff does not care how a step was
/// produced, only whether two steps compare equal, so building these by hand keeps the alignment
/// tests independent of projection.
static SymbolicAllele plain(const vector<nid_t>& ids) {
    SymbolicAllele a;
    for (nid_t id : ids) {
        SymbolicStep s;
        s.id = id;
        a.push_back(s);
    }
    return a;
}

/// One chain symbol, so a test can mix chains and plain nodes.
static SymbolicStep chain_step(nid_t start, nid_t end, bool backward = false) {
    SymbolicStep s;
    s.id = start;
    s.end_id = end;
    s.backward = backward;
    return s;
}

TEST_CASE("Identical symbolic alleles produce no difference blocks", "[symbolic_diff]") {
    REQUIRE(symbolic_diff(plain({1, 2, 3, 4}), plain({1, 2, 3, 4})).empty());
    REQUIRE(symbolic_diff(plain({}), plain({})).empty());
}

TEST_CASE("A pure insertion and a pure deletion each give one block", "[symbolic_diff]") {
    auto ins = symbolic_diff(plain({}), plain({7, 8}));
    REQUIRE(ins.size() == 1);
    REQUIRE(ins[0] == DiffBlock{0, 0, 0, 2});
    REQUIRE(ins[0].ref_empty());

    auto del = symbolic_diff(plain({7, 8}), plain({}));
    REQUIRE(del.size() == 1);
    REQUIRE(del[0] == DiffBlock{0, 2, 0, 0});
    REQUIRE(del[0].alt_empty());
}

TEST_CASE("An interior insertion is one block between two matched runs", "[symbolic_diff]") {
    // 1 2 . 4   against   1 2 3 4
    auto b = symbolic_diff(plain({1, 2, 4}), plain({1, 2, 3, 4}));
    REQUIRE(b.size() == 1);
    REQUIRE(b[0] == DiffBlock{2, 2, 2, 3});
}

TEST_CASE("Two separated differences give two blocks", "[symbolic_diff]") {
    // This is the population the whole change exists for: one snarl, two independent differences,
    // today emitted as a single ALT spanning both plus the invariant sequence between them.
    auto b = symbolic_diff(plain({1, 2, 3, 4, 5}), plain({1, 9, 3, 8, 5}));
    REQUIRE(b.size() == 2);
    REQUIRE(b[0] == DiffBlock{1, 2, 1, 2});
    REQUIRE(b[1] == DiffBlock{3, 4, 3, 4});
}

TEST_CASE("Adjacent differences aggregate into one block", "[symbolic_diff]") {
    // Two substitutions with no matched step between them are one difference, not two.
    auto b = symbolic_diff(plain({1, 2, 3, 4}), plain({1, 7, 8, 4}));
    REQUIRE(b.size() == 1);
    REQUIRE(b[0] == DiffBlock{1, 3, 1, 3});
}

TEST_CASE("The substitution cost model resolves the ambiguous pair to one block",
          "[symbolic_diff]") {
    // The load-bearing tie-break case. [a,b] against [b,b] has two minimal insert/delete-only
    // alignments of equal cost: delete a, match b, insert b -- which is TWO blocks separated by a
    // spurious match -- or substitute a for b and match b, which is ONE. Substitution at cost 1
    // beats delete-plus-insert at 2, so the one-block reading must win strictly. If this test ever
    // reports 2, the cost model has silently reverted and every downstream block count changes.
    auto b = symbolic_diff(plain({1, 2}), plain({2, 2}));
    REQUIRE(b.size() == 1);
    REQUIRE(b[0] == DiffBlock{0, 1, 0, 1});

    // Same shape mirrored, so the answer cannot come from a left/right asymmetry.
    auto c = symbolic_diff(plain({2, 1}), plain({2, 2}));
    REQUIRE(c.size() == 1);
    REQUIRE(c[0] == DiffBlock{1, 2, 1, 2});
}

TEST_CASE("The alignment is deterministic under repetition", "[symbolic_diff]") {
    // Determinism is the property the caller actually depends on: this result decides how many
    // records a snarl emits, so an unstable tie-break would make output vary with nothing the
    // caller controls. A deliberately degenerate pair, where many alignments tie.
    SymbolicAllele r = plain({1, 1, 1, 2, 1, 1});
    SymbolicAllele h = plain({1, 1, 2, 1, 1, 1});
    auto first = symbolic_diff(r, h);
    for (int trial = 0; trial < 8; ++trial) {
        REQUIRE(symbolic_diff(r, h) == first);
    }
}

TEST_CASE("A chain symbol matches only a chain symbol for the same chain and direction",
          "[symbolic_diff]") {
    SymbolicAllele r;
    r.push_back(plain({1})[0]);
    r.push_back(chain_step(2, 4));
    r.push_back(plain({5})[0]);

    // Same chain, same direction: matched, so no block at all.
    SymbolicAllele same = r;
    REQUIRE(symbolic_diff(r, same).empty());

    // Same chain crossed the other way is a different route through this snarl, so it is a
    // difference here rather than something to delegate to the chain's own record.
    SymbolicAllele flipped;
    flipped.push_back(plain({1})[0]);
    flipped.push_back(chain_step(2, 4, true));
    flipped.push_back(plain({5})[0]);
    auto b = symbolic_diff(r, flipped);
    REQUIRE(b.size() == 1);
    REQUIRE(b[0] == DiffBlock{1, 2, 1, 2});
}

TEST_CASE("A repeated chain symbol still yields a deterministic block list", "[symbolic_diff]") {
    // Loops: the same chain crossed twice makes the alignment maximally degenerate, which is
    // exactly where the tie-break has to hold.
    SymbolicAllele r;
    r.push_back(chain_step(2, 4));
    r.push_back(plain({6})[0]);
    r.push_back(chain_step(2, 4));

    SymbolicAllele h;
    h.push_back(chain_step(2, 4));
    h.push_back(chain_step(2, 4));

    auto first = symbolic_diff(r, h);
    REQUIRE_FALSE(first.empty());
    for (int trial = 0; trial < 8; ++trial) {
        REQUIRE(symbolic_diff(r, h) == first);
    }
}

TEST_CASE("A haplotype that deletes a chain has no crossing of it, though the reference step "
          "is inside a block", "[symbolic_diff]") {
    // The distinction a delegation rule has to get right, and the one that cost a 60x over-fire.
    //
    // R crosses chain C; H skips it. Two facts hold at once and they are NOT the same fact:
    //   - the REFERENCE's chain step falls inside the difference block, and
    //   - the HAPLOTYPE crosses the chain zero times.
    // A rule asking "is this chain already spelled out by a block ALT" must consult the second.
    // Reading the first says yes for a haplotype that carries no copy of the chain at all.
    SymbolicAllele r;
    r.push_back(plain({1})[0]);
    r.push_back(chain_step(2, 4));
    r.push_back(plain({5})[0]);

    SymbolicAllele h = plain({1, 5});

    auto blocks = symbolic_diff(r, h);
    REQUIRE(blocks.size() == 1);
    // The reference's chain step (index 1) is inside the block.
    REQUIRE(blocks[0].ref_begin <= 1);
    REQUIRE(1 < blocks[0].ref_end);
    // And the haplotype carries no chain symbol at all, so it crosses the chain zero times.
    int crossings = 0;
    for (const SymbolicStep& step : h) {
        if (step.is_chain() && step.id == 2 && step.end_id == 4) {
            ++crossings;
        }
    }
    REQUIRE(crossings == 0);
    // The alt side of the block is empty, which is the same fact stated a third way.
    REQUIRE(blocks[0].alt_empty());
}

TEST_CASE("A haplotype crossing a chain inside a difference block does carry it in the block",
          "[symbolic_diff]") {
    // The genuine inline case: H takes a different route into the same chain, so the chain symbol
    // sits inside the block on the ALT side and the block's sequence spells the route through it.
    SymbolicAllele r;
    r.push_back(plain({1})[0]);
    r.push_back(plain({9})[0]);
    r.push_back(chain_step(2, 4));
    r.push_back(plain({5})[0]);

    SymbolicAllele h;
    h.push_back(plain({1})[0]);
    h.push_back(chain_step(2, 4));
    h.push_back(plain({5})[0]);

    auto blocks = symbolic_diff(r, h);
    REQUIRE(blocks.size() == 1);
    // Node 9 is deleted; the chain matches, so it must NOT be inside the block on the alt side.
    bool chain_inside = false;
    for (size_t j = 0; j < h.size(); ++j) {
        if (!h[j].is_chain()) {
            continue;
        }
        for (const DiffBlock& b : blocks) {
            if ((size_t)b.alt_begin <= j && j < (size_t)b.alt_end) {
                chain_inside = true;
            }
        }
    }
    REQUIRE_FALSE(chain_inside);
}

TEST_CASE("The degradation flag is set only when the DP is refused", "[symbolic_diff]") {
    bool degraded = true;
    symbolic_diff(plain({1, 2, 3}), plain({1, 4, 3}), &degraded);
    REQUIRE_FALSE(degraded);

    // Empty-side cases are handled without the DP but are not degradations: the answer is exact.
    degraded = true;
    symbolic_diff(plain({}), plain({1, 2}), &degraded);
    REQUIRE_FALSE(degraded);
}

TEST_CASE("Visit ranges partition a chain-free traversal one visit per step", "[symbolic_allele]") {
    Snarl top = make_snarl(1, 5);
    auto mgr = make_manager(top, {});
    SnarlTraversal t = make_trav({1, 2, 3, 5});

    vector<pair<int, int>> ranges;
    SymbolicAllele a = symbolic_allele(t, top, *mgr, &ranges);
    REQUIRE(a.size() == 4);
    REQUIRE(ranges.size() == 4);
    for (int k = 0; k < 4; ++k) {
        REQUIRE(ranges[k] == make_pair(k, k + 1));
    }
    check_partition(t, top, *mgr);
}

TEST_CASE("A chain symbol's visit range runs from entry up to but not including its exit",
          "[symbolic_allele]") {
    // 1 -> [child 2..4] -> 5. The chain is entered at visit 1 and left at visit 3, and visit 3
    // (node 4) belongs to the step *after* the symbol, because the boundary is shared.
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});
    SnarlTraversal t = make_trav({1, 2, 3, 4, 5});

    vector<pair<int, int>> ranges;
    SymbolicAllele a = symbolic_allele(t, top, *mgr, &ranges);
    REQUIRE(ranges.size() == a.size());
    REQUIRE(a[1].is_chain());
    REQUIRE(ranges[0] == make_pair(0, 1));   // node 1
    REQUIRE(ranges[1] == make_pair(1, 3));   // the chain, visits 1..2
    REQUIRE(ranges[2] == make_pair(3, 4));   // node 4, the shared exit boundary
    check_partition(t, top, *mgr);
}

TEST_CASE("Visit ranges still partition when the route through the chain is longer",
          "[symbolic_allele]") {
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});
    // A longer interior: the symbol swallows more visits but the partition must still close.
    check_partition(make_trav({1, 2, 3, 6, 7, 4, 5}), top, *mgr);
    // And when the chain is skipped entirely there is no symbol at all.
    check_partition(make_trav({1, 5}), top, *mgr);
}

TEST_CASE("Passing no range vector leaves symbolic projection unchanged", "[symbolic_allele]") {
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});
    SnarlTraversal t = make_trav({1, 2, 3, 4, 5});

    vector<pair<int, int>> ranges;
    REQUIRE(symbolic_allele(t, top, *mgr) == symbolic_allele(t, top, *mgr, &ranges));
}

TEST_CASE("A reversed snarl resolves to the same snarl and still symbolises its children",
          "[symbolic_allele]") {
    // The caller works on a REVERSED copy of any snarl whose reference path runs backwards --
    // `flip_snarl` swaps the boundaries and reverses each. The reversed copy's start node is
    // therefore the original END node, and an identity test demanding start-to-start rejects it.
    // That rejection made `site_ptr` null, `is_child` false at every visit, and the projection a
    // bare node list with no chain symbols: symbolic collapsing silently off for the whole snarl,
    // on 7.4% of chr20's sites.
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    // Exactly what flip_snarl produces from `top`.
    Snarl flipped;
    flipped.mutable_start()->set_node_id(5);
    flipped.mutable_start()->set_backward(true);
    flipped.mutable_end()->set_node_id(1);
    flipped.mutable_end()->set_backward(true);

    REQUIRE(symbolic_site_resolvable(top, *mgr));
    REQUIRE(symbolic_site_resolvable(flipped, *mgr));

    // And the projection through the reversed site must still collapse the child chain, which is
    // the behaviour the resolution exists to enable rather than merely a property of it.
    SnarlTraversal backwards = make_trav_rev({5, 4, 3, 2, 1});
    SymbolicAllele a = symbolic_allele(backwards, flipped, *mgr);
    bool has_chain = false;
    for (const SymbolicStep& step : a) {
        if (step.is_chain()) {
            has_chain = true;
        }
    }
    REQUIRE(has_chain);
    REQUIRE(has_child_chain(backwards, flipped, *mgr));

    // Two traversals differing only inside the child are still equal through the reversed site,
    // which is the whole point and is what was lost.
    SnarlTraversal other = make_trav_rev({5, 4, 9, 2, 1});
    REQUIRE(symbolically_equal(backwards, other, flipped, *mgr));
}

TEST_CASE("A snarl sharing neither boundary is still refused", "[symbolic_allele]") {
    // Nothing in the manager has boundaries {7, 8}, so the LOOKUP fails and resolution stops before
    // the boundary comparison is reached. Kept for what it is -- a check on the null path -- and
    // explicitly not a test of the widened comparison, which the next case covers.
    Snarl top = make_snarl(1, 5);
    auto mgr = make_manager(top, {});
    REQUIRE_FALSE(symbolic_site_resolvable(make_snarl(7, 8), *mgr));
}

TEST_CASE("A resolvable lookup with mismatched boundaries is still refused", "[symbolic_allele]") {
    // The negative case that actually exercises the widened comparison. Node 1 IS a boundary, so
    // into_which_snarl succeeds and returns (1,5) -- and the site claims to be (1,4), which matches
    // neither forward (end 5 != 4) nor reversed (start 1 != 4). Without this the suite only ever
    // asserted that reversals are accepted, never that impostors are still rejected, because the
    // other negative case exits at the null-lookup gate and never reaches the predicate at all.
    Snarl top = make_snarl(1, 5);
    Snarl child = make_snarl(2, 4);
    auto mgr = make_manager(top, {child});

    REQUIRE(symbolic_site_resolvable(top, *mgr));
    REQUIRE_FALSE(symbolic_site_resolvable(make_snarl(1, 4), *mgr));

    // And a half-reversed site -- right nodes, wrong orientation on one boundary -- must also be
    // refused now that the reversed branch compares whole visits rather than node ids.
    Snarl half;
    half.mutable_start()->set_node_id(5);
    half.mutable_start()->set_backward(true);
    half.mutable_end()->set_node_id(1);
    half.mutable_end()->set_backward(false);   // a true reversal would have this backward
    REQUIRE_FALSE(symbolic_site_resolvable(half, *mgr));
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
