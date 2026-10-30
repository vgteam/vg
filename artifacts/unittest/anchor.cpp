/// \file unittest/anchor.cpp
///
/// Unit tests for anchor pin resolution.
///
/// The pin is a zero-length point between two positions, so everything here is off-by-one and
/// strand arithmetic. These tests pin the convention; the check made while pins are resolved
/// (`AnchorPinCounters::verify_failed`) catches errors on real data.
///

#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../anchor.hpp"
#include "../../src/unittest/catch.hpp"
#include "../../src/snarls.hpp"
#include "../../src/utility.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A site with a SNP and a deletion of the SNP-bearing node:
///
///        2 (T)
///      /      \
///   1 --- 3 (G) --- 4          and the deletion edge 1 -> 4
///      \__________/
///
/// Node 1 is "AAAACCCC" and node 4 "GGGGTTTT", both 8 bp, so a pin at either boundary has a
/// distinctive base on each side.
struct PinSite {
    bdsg::HashGraph graph;
    Snarl snarl;

    PinSite() {
        graph.create_handle("AAAACCCC", 1);
        graph.create_handle("T", 2);
        graph.create_handle("G", 3);
        graph.create_handle("GGGGTTTT", 4);
        graph.create_edge(graph.get_handle(1), graph.get_handle(2));
        graph.create_edge(graph.get_handle(2), graph.get_handle(4));
        graph.create_edge(graph.get_handle(1), graph.get_handle(3));
        graph.create_edge(graph.get_handle(3), graph.get_handle(4));
        graph.create_edge(graph.get_handle(1), graph.get_handle(4));
        snarl.mutable_start()->set_node_id(1);
        snarl.mutable_end()->set_node_id(4);
        snarl.set_type(ULTRABUBBLE);
        snarl.set_start_end_reachable(true);
    }
};

/// An all-match alignment over the given (node, is_reverse) steps.
static Alignment matching_read(const HandleGraph& graph, const string& name,
                               const vector<pair<nid_t, bool>>& steps) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (const auto& step : steps) {
        string node_seq = graph.get_sequence(graph.get_handle(step.first, step.second));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(step.first);
        m->mutable_position()->set_is_reverse(step.second);
        m->mutable_position()->set_offset(0);
        Edit* e = m->add_edit();
        e->set_from_length(node_seq.size());
        e->set_to_length(node_seq.size());
        seq += node_seq;
    }
    aln.set_sequence(seq);
    aln.set_mapping_quality(60);
    return aln;
}

TEST_CASE("A forward read pins at the base before each junction", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // 1 -> 2 -> 4, all forward: "AAAACCCC" "T" "GGGGTTTT".
    Alignment aln = matching_read(site.graph, "fwd", {{1, false}, {2, false}, {4, false}});

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(s.placed());
    // The pin follows node 1's last base, which is read index 7.
    REQUIRE(s.offset == 7);
    REQUIRE(s.direction == 0);

    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(e.placed());
    // The pin precedes node 4's first base at read index 9, so the last base BEFORE it is node 2's
    // single base at index 8 -- inside the site, which is where the entry pin's upstream side is.
    REQUIRE(e.offset == 8);
    REQUIRE(e.direction == 0);
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("A reverse read pins on the same graph positions with strand 1", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // The same walk read the other way: 4rev -> 2rev -> 1rev.
    Alignment aln = matching_read(site.graph, "rev", {{4, true}, {2, true}, {1, true}});
    REQUIRE(aln.sequence() == "AAAACCCCAGGGGTTTT");

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(s.placed());
    // Node 1's last base in the SITE's direction is its node-forward index 7, which this read places
    // at index 9 -- the first base of its node-1 mapping, because it reads the node backwards.
    REQUIRE(s.offset == 9);
    REQUIRE(s.direction == 1);
    // And the read base there is the complement of the graph's, which is what the invariant checks.
    REQUIRE(aln.sequence()[9] == 'G');

    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(e.placed());
    // Node 4's first base in the site's direction sits at read index 7; upstream in the site's
    // direction is LATER in this read, so the offset steps forward rather than back.
    REQUIRE(e.offset == 8);
    REQUIRE(e.direction == 1);
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("A soft clip shifts the offset without moving the pin", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    Alignment aln = matching_read(site.graph, "clipped", {{1, false}, {2, false}, {4, false}});
    // Three unaligned bases at the front: a soft clip is an edit that consumes read but no node.
    Mapping* first = aln.mutable_path()->mutable_mapping(0);
    Edit clip;
    clip.set_from_length(0);
    clip.set_to_length(3);
    clip.set_sequence("GGG");
    *first->add_edit() = *first->mutable_edit(0);
    *first->mutable_edit(0) = clip;
    aln.set_sequence("GGG" + aln.sequence());

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(s.placed());
    REQUIRE(s.offset == 10);
    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(e.placed());
    REQUIRE(e.offset == 11);
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("An insertion at the junction falls after the pin", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    Alignment aln = matching_read(site.graph, "ins", {{1, false}, {2, false}, {4, false}});
    // Two inserted bases after node 2's base, so they sit between the site interior and node 4.
    Mapping* mid = aln.mutable_path()->mutable_mapping(1);
    Edit* ins = mid->add_edit();
    ins->set_from_length(0);
    ins->set_to_length(2);
    ins->set_sequence("CC");
    aln.set_sequence("AAAACCCCT" "CC" "GGGGTTTT");

    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(e.placed());
    // Node 4's first base is now at read index 11, but indices 9 and 10 consume no node base, so
    // the pin steps over them to node 2's base at 8 rather than reporting an inserted base.
    REQUIRE(e.offset == 8);
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("A mismatch at the pin is still placeable and still verified", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    Alignment aln = matching_read(site.graph, "snv", {{1, false}, {2, false}, {4, false}});
    // Node 1's last base wrong: split its 8M into 7M + 1X.
    Mapping* first = aln.mutable_path()->mutable_mapping(0);
    first->clear_edit();
    Edit* match = first->add_edit();
    match->set_from_length(7);
    match->set_to_length(7);
    Edit* sub = first->add_edit();
    sub->set_from_length(1);
    sub->set_to_length(1);
    sub->set_sequence("A");
    aln.set_sequence("AAAACCCA" "T" "GGGGTTTT");

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(s.placed());
    REQUIRE(s.offset == 7);
    // The base legitimately differs from the graph, so the edit's own sequence is the reference the
    // invariant uses instead. It still checks the offset arithmetic.
    REQUIRE(counters.pins.verified.load() == 1);
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("A deletion over the pin's node base refuses the pin rather than walking back",
          "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    Alignment aln = matching_read(site.graph, "del", {{1, false}, {2, false}, {4, false}});
    // Node 1's last base deleted: 7M then 1D.
    Mapping* first = aln.mutable_path()->mutable_mapping(0);
    first->clear_edit();
    Edit* match = first->add_edit();
    match->set_from_length(7);
    match->set_to_length(7);
    Edit* del = first->add_edit();
    del->set_from_length(1);
    del->set_to_length(0);
    aln.set_sequence("AAAACCC" "T" "GGGGTTTT");

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    // Refused, NOT resolved to the previous base. On a 1 bp boundary node walking back would step
    // straight past the node and land on the neighbouring snarl's pin position.
    REQUIRE(!s.placed());
    REQUIRE(counters.pins.unaligned_base.load() == 1);
}

TEST_CASE("A read that does not reach the junction is refused", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // Starts at node 4, so nothing of it lies upstream of the entry pin.
    Alignment aln = matching_read(site.graph, "starts_here", {{4, false}});
    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(!e.placed());
    REQUIRE(counters.pins.no_neighbour.load() == 1);

    // And a read that never visits the node at all.
    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(!s.placed());
    REQUIRE(counters.pins.no_visit.load() == 1);
}

TEST_CASE("A read spanning the deletion edge lands both pins on one position", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // 1 -> 4 directly: the allele that deletes the whole site, so the read has no bases inside it.
    Alignment aln = matching_read(site.graph, "spanner", {{1, false}, {4, false}});

    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(s.placed());
    REQUIRE(e.placed());
    // This is the one case where two pins of one snarl coincide in a read, and it is why
    // build_site_anchors keeps such a read at S and drops it from E.
    REQUIRE(s.offset == e.offset);
    REQUIRE(s.direction == e.direction);
}

TEST_CASE("A deleted node upstream of an entry pin refuses it, rather than stepping over it",
          "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // 1 -> 2 -> 4, but node 2's single base is deleted in this read. The read still WALKS through
    // node 2 -- it just places no base there.
    Alignment aln = matching_read(site.graph, "gap", {{1, false}, {2, false}, {4, false}});
    Mapping* mid = aln.mutable_path()->mutable_mapping(1);
    mid->clear_edit();
    Edit* del = mid->add_edit();
    del->set_from_length(1);
    del->set_to_length(0);
    aln.set_sequence("AAAACCCC" "GGGGTTTT");

    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    // Node 4's first base is at read index 8, and the base before it in READ order is node 1's last
    // base at index 7 -- but that is a whole node further upstream in the GRAPH, and index 7 is
    // where the snarl's own start pin sits. Reporting it would put one read position in two anchors,
    // which is the thing the pin geometry exists to prevent. So the pin is refused.
    REQUIRE(!e.placed());
    REQUIRE(counters.pins.no_neighbour.load() == 1);

    // The start pin is unaffected: its reference base is inside its own node.
    AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters);
    REQUIRE(s.placed());
    REQUIRE(s.offset == 7);
}

TEST_CASE("An entry pin steps one node, not one aligned base", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // Node 2 present but carrying an insertion after its base, so the read base immediately before
    // node 4 is inserted. The pin must report node 2's own base, not the inserted one.
    Alignment aln = matching_read(site.graph, "ins2", {{1, false}, {2, false}, {4, false}});
    Mapping* mid = aln.mutable_path()->mutable_mapping(1);
    Edit* ins = mid->add_edit();
    ins->set_from_length(0);
    ins->set_to_length(3);
    ins->set_sequence("AAA");
    aln.set_sequence("AAAACCCCT" "AAA" "GGGGTTTT");

    AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 4, false, false, counters);
    REQUIRE(e.placed());
    REQUIRE(e.offset == 8);        // node 2's base, not index 11
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

/// Evidence for one site with two alleles, `n` reads, each perfectly fitting allele `prefers[i]`.
static AnchorSiteEvidence two_allele_evidence(const vector<int>& prefers,
                                              const vector<AnchorPlacement>& start_pins,
                                              const vector<AnchorPlacement>& end_pins) {
    AnchorSiteEvidence ev;
    ev.n_alleles = 2;
    ev.start_node = 1;
    ev.end_node = 4;
    ev.allele_length = {9, 9};
    ev.mean_read_length = 10.0f;
    ev.length_weighted = true;
    for (size_t i = 0; i < prefers.size(); ++i) {
        AnchorRead r;
        r.read = read_names().intern("r" + std::to_string(i));
        r.mismap = 0.02f;
        r.start_pin = start_pins[i];
        r.end_pin = end_pins[i];
        ev.reads.push_back(r);
        ev.rel.push_back(prefers[i] == 0 ? 1.0f : 0.0f);
        ev.rel.push_back(prefers[i] == 1 ? 1.0f : 0.0f);
    }
    return ev;
}

TEST_CASE("A coincident read is kept at one pin only", "[anchor]") {
    // Four reads, two per allele; every one of them has both pins at the same position, which is
    // what a read spanning a whole-site deletion looks like.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement p;
        p.offset = 10 + i;
        p.direction = 0;
        starts.push_back(p);
        ends.push_back(p);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 1, 1}, starts, ends);

    AnchorParams params;
    params.min_reads = 1;
    AnchorCounters counters;
    params.counters = &counters;
    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);

    REQUIRE(counters.coincident.load() == 4);
    // Only the two start anchors survive; every end anchor lost all its reads.
    REQUIRE(out.size() == 2);
    for (const auto& anchor : out) {
        REQUIRE(anchor.node == 1);
    }
}

TEST_CASE("The two slots of a site partition the reads, and score is mismap-bounded", "[anchor]") {
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        s.offset = 10 + i;
        starts.push_back(s);
        AnchorPlacement e;
        e.offset = 30 + i;   // distinct from the start pin, so both survive
        ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 1, 1}, starts, ends);

    AnchorParams params;
    params.min_reads = 1;
    AnchorCounters counters;
    params.counters = &counters;
    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);

    REQUIRE(out.size() == 4);   // two slots at each of two pins
    size_t total = 0;
    for (const auto& anchor : out) {
        REQUIRE(anchor.reads.size() == 2);
        total += anchor.reads.size();
        for (const auto& row : anchor.reads) {
            // A perfectly discriminating read on a balanced het, at a 0.02 mismap floor, is worth
            // about 14 phred and no more. The bound is the point: it is what stops a per-read score
            // claiming more confidence than the read's own mapping supports.
            REQUIRE(row.score > 13.0);
            REQUIRE(row.score < 15.0);
        }
    }
    REQUIRE(total == 8);
}

TEST_CASE("An end pin holding no reads of its own is dropped", "[anchor]") {
    // Four reads, every one of them placed at BOTH pins. The end pin therefore holds nothing the
    // start pin does not, and the two are joined by all the same reads.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s; s.offset = 10 + i; starts.push_back(s);
        AnchorPlacement e; e.offset = 30 + i; ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 1, 1}, starts, ends);
    AnchorParams params;
    params.min_reads = 1;
    AnchorCounters counters;
    params.counters = &counters;

    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 4);            // two slots at each of two pins

    out.clear();
    params.end_pin_min_new = 1;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 2);            // two slots, start pin only
    for (const auto& anchor : out) {
        REQUIRE(anchor.node == ev.start_node);
    }
    // Counted per slot, because the decision is per slot: both end anchors were redundant.
    REQUIRE(counters.single_pin.load() == 2);
}

TEST_CASE("An end pin reached by a read the start pin misses is kept", "[anchor]") {
    // One read reaches only the end pin -- it started inside the site, so it never crossed the start
    // boundary. That read's linkage exists nowhere else, so the end pin has to stay.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        if (i != 3) {
            s.offset = 10 + i;           // read 3 does not reach the start pin
        }
        starts.push_back(s);
        AnchorPlacement e; e.offset = 30 + i; ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 1, 1}, starts, ends);
    AnchorParams params;
    params.min_reads = 1;
    params.end_pin_min_new = 1;
    AnchorCounters counters;
    params.counters = &counters;

    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);

    // Read 3 sits in slot 1, so it is slot 1's end anchor that earns its place. Slot 0's end anchor
    // is entirely redundant and goes -- the decision is per slot, and a site-level one would have
    // kept a wholly redundant anchor on the strength of a read in the other slot.
    REQUIRE(out.size() == 3);
    REQUIRE(counters.single_pin.load() == 1);
    size_t end_anchors = 0;
    for (const auto& anchor : out) {
        if (anchor.node == ev.end_node) {
            ++end_anchors;
            REQUIRE(anchor.slot == 1);
        }
    }
    REQUIRE(end_anchors == 1);

    // Raising the bar past what it contributes drops it too: one new read is not always worth
    // carrying a whole anchor's worth of otherwise-redundant rows.
    out.clear();
    params.end_pin_min_new = 2;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 2);
}

TEST_CASE("A homozygous site is emitted unless excluded", "[anchor]") {
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        s.offset = 10 + i;
        starts.push_back(s);
        AnchorPlacement e;
        e.offset = 30 + i;
        ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 0, 0}, starts, ends);

    AnchorCounters counters;
    AnchorParams params;
    params.counters = &counters;
    params.min_reads = 1;
    vector<AnchorWriter::Anchor> out;
    // By default a homozygous site still anchors: one anchor per pin, both slots collapsed into one.
    // It partitions nothing, but it links reads, and an anchor graph needs contiguity too.
    build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 2);
    REQUIRE(out[0].slot == 0);
    REQUIRE(out[0].reads.size() == 4);

    out.clear();
    params.het_only = true;
    build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.empty());
}

TEST_CASE("A homozygous site splits by read phase only when both strands are supported",
          "[anchor]") {
    // Both haplotypes at a homozygous site spell the same allele, so the site's alleles say
    // nothing about which strand a read is on. With `hom_split` on, the site's reads are divided
    // between two slots by the strand log-odds they carry from other sites. The site is split
    // only when its reads lean both ways. If every read leaned one way, they would all land in
    // one slot, which relabels the site without dividing its reads.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        s.offset = 10 + i;
        starts.push_back(s);
        AnchorPlacement e;
        e.offset = 30 + i;
        ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 0, 0}, starts, ends);

    AnchorCounters counters;
    AnchorParams params;
    params.counters = &counters;
    params.min_reads = 1;
    params.hom_split = true;
    // Pinned, not inherited. The fixture has four reads, so the shipped default could never be met
    // and every section below would collapse for the wrong reason -- the two that expect a
    // collapse would stop testing one-sidedness and start testing an unreachable quota, which is
    // the shape of a check that cannot fire.
    params.phase_min_side = 2;
    vector<AnchorWriter::Anchor> out;

    SECTION("reads that lean both ways are split into two slots carrying the SAME allele") {
        vector<double> strand{+5.0, +5.0, -5.0, -5.0};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        // Two pins x two slots.
        REQUIRE(out.size() == 4);
        size_t in_slot[2] = {0, 0};
        for (const AnchorWriter::Anchor& a : out) {
            REQUIRE(a.slot >= 0);
            REQUIRE(a.slot <= 1);
            REQUIRE(a.allele == 0);   // equal alleles across the two slots is the tell
            in_slot[a.slot] += a.reads.size();
        }
        // Each pin carries the same partition, so 2 reads per slot per pin.
        REQUIRE(in_slot[0] == 4);
        REQUIRE(in_slot[1] == 4);
    }

    SECTION("reads that all lean one way leave the site collapsed") {
        vector<double> strand{+5.0, +5.0, +5.0, +5.0};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(out.size() == 2);       // one anchor per pin, as with `hom_split` off
        REQUIRE(out[0].slot == 0);
        REQUIRE(out[0].reads.size() == 4);
    }

    SECTION("reads with no cross-site opinion leave the site collapsed") {
        vector<double> strand{0.0, 0.0, 0.0, 0.0};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(out.size() == 2);
        REQUIRE(out[0].reads.size() == 4);
    }

    SECTION("a read with no cross-site opinion is placed by the coin, on ONE haplotype") {
        // `lo > 0 ? 0 : 1` would send every zero to slot 1 -- systematically, in one direction,
        // with nothing behind it. Dropping it instead was the first answer and it was wrong for a
        // different reason: the same read IS placed at het sites and at unsplit homozygous ones, so
        // leaving it out here alone punches a hole in its trail, and a read that vanishes between
        // two anchors and comes back is the shape of a deletion.
        //
        // So it is placed, by the read name's hash parity. The site is homozygous, so both slots
        // spell the same allele and the read's sequence fits either haplotype: only the label is
        // arbitrary. What must hold is that the label is CONSISTENT -- one haplotype, at both pins.
        // The earlier all-zero section cannot catch any of this: there the site never splits.
        params.phase_min_side = 1;   // so three opinionated reads are enough to split
        vector<double> strand{+5.0, -5.0, +5.0, 0.0};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(out.size() == 4);
        size_t placed = 0;
        for (const AnchorWriter::Anchor& a : out) {
            placed += a.reads.size();
        }
        // All 4 reads x 2 pins: no read is dropped for want of strand log-odds.
        REQUIRE(placed == 8);
        REQUIRE(counters.hom_split_no_opinion.load() > 0);
        REQUIRE(counters.hom_split_coin.load() > 0);

        // The coin's load-bearing property: the no-opinion read is on exactly one slot, and the
        // same one at both pins. A per-site or per-pin random choice would satisfy the count above
        // and fail this.
        const uint32_t quiet = ev.reads[3].read;
        set<int> slots_held;
        size_t rows = 0;
        for (const AnchorWriter::Anchor& a : out) {
            for (const AnchorWriter::ReadRow& row : a.reads) {
                if (row.read == quiet) {
                    slots_held.insert(a.slot);
                    ++rows;
                }
            }
        }
        REQUIRE(rows == 2);                 // one row per pin
        REQUIRE(slots_held.size() == 1);    // and both on the same haplotype
    }

    SECTION("a read whose strand is from another phase set is left out of a split site only") {
        // NaN is what `read_strand_log_odds` returns for a read whose strand belongs to another
        // phase set, or to more than one. Its strand says nothing about this phase set's, so a
        // split site drops it rather than letting the coin pick a strand for it; where it also
        // crosses a heterozygous site of this phase set, a coin would contradict the slot that
        // site's allele gives it half the time.
        params.phase_min_side = 1;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        vector<double> strand{+5.0, -5.0, +5.0, nan};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(out.size() == 4);
        size_t placed = 0;
        for (const AnchorWriter::Anchor& a : out) {
            placed += a.reads.size();
            for (const AnchorWriter::ReadRow& row : a.reads) {
                REQUIRE(row.read != ev.reads[3].read);
            }
        }
        REQUIRE(placed == 6);   // three reads x two pins
        REQUIRE(counters.hom_split_no_opinion.load() == 1);
        REQUIRE(counters.hom_split_coin.load() == 0);

        // At a site left collapsed no strand is used, so the same read is kept there.
        out.clear();
        vector<double> one_sided{+5.0, +5.0, +5.0, nan};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &one_sided);
        REQUIRE(out.size() == 2);
        REQUIRE(out[0].reads.size() == 4);
    }

    SECTION("with hom_split off the same evidence stays collapsed") {
        vector<double> strand{+5.0, +5.0, -5.0, -5.0};
        params.hom_split = false;
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(out.size() == 2);
        REQUIRE(out[0].reads.size() == 4);
    }

    SECTION("a split read keeps the score it would have had unsplit") {
        // Both slots spell one allele, so the share must be taken over that one allele. Halving
        // the mixture would drop a clean read from about 13 phred to about 3. That would change
        // what `min_read_score` and the reliability column mean at split sites.
        vector<AnchorWriter::Anchor> unsplit;
        params.hom_split = false;
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, unsplit);
        REQUIRE(!unsplit.empty());
        const float want = unsplit[0].reads.at(0).score;

        out.clear();
        params.hom_split = true;
        vector<double> strand{+5.0, +5.0, -5.0, -5.0};
        build_site_anchors(ev, AnchorCall{{0, 0}, 0, 0.9, 1.0}, ">1>4", params, out, &strand);
        REQUIRE(!out.empty());
        REQUIRE(out[0].reads.at(0).score == Approx(want));
    }
}

TEST_CASE("A site's reliability is the mean score of the reads it emitted", "[anchor]") {
    // Site-level, so every anchor of the site carries the same value, and it is computed over the
    // reads that actually reached the file -- a consumer averaging the R rows must get the same
    // number back, which is the only reason to write a derivable column at all.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        s.offset = 10 + i;
        starts.push_back(s);
        AnchorPlacement e;
        e.offset = 30 + i;
        ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 1, 1}, starts, ends);

    AnchorCounters counters;
    AnchorParams params;
    params.counters = &counters;
    params.min_reads = 1;
    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 4);

    // Every row of the site agrees, and none is the "no reads" sentinel.
    for (const AnchorWriter::Anchor& a : out) {
        REQUIRE(a.reliability == Approx(out[0].reliability));
        REQUIRE(a.reliability >= 0.0);
    }

    // And it equals the mean over the site's DISTINCT reads -- each counted once, though every
    // one of them appears at two pins.
    std::map<string, double> per_read;
    size_t rows = 0;
    for (const AnchorWriter::Anchor& a : out) {
        for (const AnchorWriter::ReadRow& r : a.reads) {
            per_read[string(read_names().name(r.read))] = r.score;
            ++rows;
        }
    }
    REQUIRE(per_read.size() == 4);
    REQUIRE(rows == 8);                  // four reads, two pins each: the dedupe is not a no-op
    double sum = 0.0;
    for (const auto& e : per_read) {
        sum += e.second;
    }
    REQUIRE(out[0].reliability == Approx(sum / 4.0));
}

TEST_CASE("A nested haploid site takes the slot its strand names", "[anchor]") {
    // One allele, not two: a nested chain the parent's other allele deletes, so there is nothing to
    // genotype on the other strand. It has one slot like a homozygote, but that slot is a strand --
    // the VCF writes the site as `a|.` or `.|a` from the same strand -- so it must not always be
    // slot 0.
    vector<AnchorPlacement> starts, ends;
    for (int i = 0; i < 4; ++i) {
        AnchorPlacement s;
        s.offset = 10 + i;
        starts.push_back(s);
        AnchorPlacement e;
        e.offset = 30 + i;
        ends.push_back(e);
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 0, 0, 0}, starts, ends);

    AnchorCounters counters;
    AnchorParams params;
    params.counters = &counters;
    params.min_reads = 1;

    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0}, 0, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 2);            // one per pin
    REQUIRE(out[0].slot == 0);
    REQUIRE(out[1].slot == 0);
    REQUIRE(out[0].reads.size() == 4);   // the single slot holds every read either way

    out.clear();
    build_site_anchors(ev, AnchorCall{{0}, 1, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 2);
    REQUIRE(out[0].slot == 1);
    REQUIRE(out[1].slot == 1);
    REQUIRE(out[0].reads.size() == 4);   // strand changes the label, not the partition

    // The strand is meaningless for a pair and must not leak into one: a diploid site's slots are
    // its GT's field order, which `phase_ordered_genotype` has already applied. Its own evidence,
    // with the reads actually split between the alleles -- against the all-allele-0 fixture above,
    // slot 1 would hold nothing and `min_reads` would drop it before the slot could be checked.
    AnchorSiteEvidence split = two_allele_evidence({0, 0, 1, 1}, starts, ends);
    out.clear();
    build_site_anchors(split, AnchorCall{{0, 1}, 1, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.size() == 4);
    int seen0 = 0, seen1 = 0;
    for (const AnchorWriter::Anchor& a : out) {
        REQUIRE(a.slot >= 0);
        REQUIRE(a.slot <= 1);
        a.slot == 0 ? ++seen0 : ++seen1;
    }
    REQUIRE(seen0 == 2);
    REQUIRE(seen1 == 2);

    // With `het_only` on, the haploid site is dropped whatever its strand, since `het_only` keeps
    // only heterozygous sites. The strand changes only which slot holds the reads.
    out.clear();
    params.het_only = true;
    build_site_anchors(ev, AnchorCall{{0}, 1, 0.9, 1.0}, ">1>4", params, out);
    REQUIRE(out.empty());
}

TEST_CASE("A pin inside a node is resolved against the base before it or the base after it",
          "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    // Pin 5 of node 1 "AAAACCCC" lies between its bases 4 and 5.
    Alignment aln = matching_read(site.graph, "fwd", {{1, false}, {2, false}, {4, false}});

    SECTION("an exit pin against base 4 and an entry pin against base 5 meet at one point") {
        AnchorPlacement s = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, true, counters, 5);
        AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, false, counters, 5);
        REQUIRE(s.offset == 4);
        REQUIRE(e.offset == 4);
        REQUIRE(e.direction == 0);
    }
    SECTION("a reverse read crosses the entry pin on strand 1, after base 5 in the read") {
        Alignment rev = matching_read(site.graph, "rev", {{4, true}, {2, true}, {1, true}});
        // Node 1 reversed starts at read index 9, so its base b is at index 9 + 7 - b.
        AnchorPlacement e = resolve_anchor_pin(SiteRead{&rev}, site.graph, 1, false, false, counters, 5);
        REQUIRE(e.direction == 1);
        REQUIRE(e.offset == 9 + 7 - 5 + 1);
    }
    SECTION("bases inserted before base 5 lie before an entry pin's point") {
        Mapping* first = aln.mutable_path()->mutable_mapping(0);
        first->clear_edit();
        for (auto lengths : vector<pair<int, int>>{{5, 5}, {0, 2}, {3, 3}}) {
            Edit* e = first->add_edit();
            e->set_from_length(lengths.first);
            e->set_to_length(lengths.second);
            if (lengths.first == 0) {
                e->set_sequence("GG");
            }
        }
        aln.set_sequence("AAAAC" "GG" "CCC" "T" "GGGGTTTT");
        AnchorPlacement e = resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, false, counters, 5);
        REQUIRE(e.offset == 6);
    }
    SECTION("a read deleting base 5 is not placed at an entry pin before it") {
        Mapping* first = aln.mutable_path()->mutable_mapping(0);
        first->mutable_edit(0)->set_from_length(5);
        first->mutable_edit(0)->set_to_length(5);
        Edit* e = first->add_edit();
        e->set_from_length(1);
        e = first->add_edit();
        e->set_from_length(2);
        e->set_to_length(2);
        aln.set_sequence("AAAAC" "CC" "T" "GGGGTTTT");
        REQUIRE(!resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, false, counters, 5).placed());
        REQUIRE(counters.pins.unaligned_base.load() == 1);
    }
    SECTION("a read starting at base 5 has no base before an entry pin there") {
        Mapping* first = aln.mutable_path()->mutable_mapping(0);
        first->mutable_position()->set_offset(5);
        first->mutable_edit(0)->set_from_length(3);
        first->mutable_edit(0)->set_to_length(3);
        aln.set_sequence("CCC" "T" "GGGGTTTT");
        REQUIRE(!resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, false, counters, 5).placed());
        REQUIRE(counters.pins.no_neighbour.load() == 1);
    }
    SECTION("a read soft clipped before base 5 has no aligned base before an entry pin there") {
        Mapping* first = aln.mutable_path()->mutable_mapping(0);
        first->mutable_position()->set_offset(5);
        first->clear_edit();
        Edit* clip = first->add_edit();
        clip->set_to_length(2);
        clip->set_sequence("GG");
        Edit* match = first->add_edit();
        match->set_from_length(3);
        match->set_to_length(3);
        aln.set_sequence("GG" "CCC" "T" "GGGGTTTT");
        REQUIRE(!resolve_anchor_pin(SiteRead{&aln}, site.graph, 1, false, false, counters, 5).placed());
        REQUIRE(counters.pins.no_neighbour.load() == 1);
    }
    REQUIRE(counters.pins.verify_failed.load() == 0);
}

TEST_CASE("A placement names the alignment, which mates sharing a name do not share", "[anchor]") {
    PinSite site;
    AnchorCounters counters;
    Alignment one = matching_read(site.graph, "pair", {{1, false}, {2, false}, {4, false}});
    Alignment other = matching_read(site.graph, "pair", {{4, true}, {2, true}, {1, true}});
    const AnchorPlacement a = resolve_anchor_pin(SiteRead{&one}, site.graph, 1, false, true, counters);
    const AnchorPlacement b = resolve_anchor_pin(SiteRead{&other}, site.graph, 1, false, true, counters);
    REQUIRE(a.alignment != b.alignment);
    REQUIRE(a.alignment == resolve_anchor_pin(SiteRead{&one}, site.graph, 4, false, false, counters).alignment);
}

TEST_CASE("An off-panel site whose two pins are one point writes one anchor per slot", "[anchor]") {
    vector<AnchorPlacement> starts(4), ends(4);
    for (int i = 0; i < 4; ++i) {
        starts[i].offset = 10 + i;
        ends[i].offset = 20 + i;
    }
    AnchorSiteEvidence ev = two_allele_evidence({0, 1, 0, 1}, starts, ends);
    ev.inside_node = true;
    ev.start_node = 7;
    ev.end_node = 7;
    ev.start_pin = 3;
    ev.end_pin = 3;
    AnchorParams params;
    params.min_reads = 1;
    AnchorCounters counters;
    params.counters = &counters;
    vector<AnchorWriter::Anchor> out;
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, "op7.2.IG", params, out);
    REQUIRE(out.size() == 2);
    REQUIRE(counters.degenerate_site.load() == 1);
    for (const AnchorWriter::Anchor& a : out) {
        REQUIRE(a.edit);
        REQUIRE(a.pin == 3);
    }
    // Pins at two points of one node are two anchors per slot.
    ev.end_pin = 4;
    out.clear();
    build_site_anchors(ev, AnchorCall{{0, 1}, 0, 0.9, 1.0}, "op7.3.SG", params, out);
    REQUIRE(out.size() == 4);
}

/// Write `anchors` with `min_reads`, and return each A row's site with the number of read rows under
/// it, in file order.
static vector<string> write_anchors(vector<AnchorWriter::Anchor> anchors, size_t min_reads,
                                    AnchorWriter& writer) {
    for (AnchorWriter::Anchor& a : anchors) {
        writer.add(std::move(a));
    }
    const string path = temp_file::create();
    AnchorParams params;
    params.min_reads = min_reads;
    REQUIRE(writer.write(path, "g", "s", "r", 0.01, params));
    vector<pair<string, size_t>> rows;
    ifstream in(path);
    string line;
    while (getline(in, line)) {
        if (line[0] == 'A') {
            istringstream fields(line);
            string kind, node, snarl;
            fields >> kind >> node >> snarl;
            rows.emplace_back(snarl, 0);
        } else if (line[0] == 'R') {
            ++rows.back().second;
        }
    }
    temp_file::remove(path);
    vector<string> out;
    for (const auto& row : rows) {
        out.push_back(row.first + ":" + to_string(row.second));
    }
    return out;
}

/// An anchor at `node` of site `snarl`, with a read row for each (read, alignment, direction,
/// offset).
static AnchorWriter::Anchor anchor_of(nid_t node, const string& snarl, bool edit,
                                      const vector<tuple<int, uint32_t, uint8_t, int64_t>>& rows) {
    AnchorWriter::Anchor a;
    a.node = node;
    a.snarl = snarl;
    a.edit = edit;
    for (const auto& r : rows) {
        AnchorWriter::ReadRow row;
        row.read = read_names().intern("point" + to_string(get<0>(r)));
        row.alignment = get<1>(r);
        row.direction = get<2>(r);
        row.offset = get<3>(r);
        a.reads.push_back(row);
    }
    return a;
}

TEST_CASE("The writer keeps each read point in one anchor, a snarl's where it has one", "[anchor]") {
    AnchorWriter writer(1);
    // Strand 0 offset 9 and strand 1 offset 10 are both the point before base 10 of a read.
    SECTION("a snarl's anchor keeps a point an off-panel anchor written before it holds") {
        const auto rows = write_anchors({anchor_of(2, "op2.4.SA", true, {{0, 1, 1, 10}, {1, 1, 0, 50}}),
                                         anchor_of(5, ">5>6", false, {{0, 1, 0, 9}})},
                                        1, writer);
        REQUIRE(rows == vector<string>{"op2.4.SA:1", ">5>6:1"});
        REQUIRE(writer.tally().placements_left_out == 1);
        REQUIRE(writer.tally().placements_written == 2);
    }
    SECTION("of two off-panel anchors, the first in the file keeps it") {
        const auto rows = write_anchors({anchor_of(3, "op3.0.SA", true, {{0, 1, 1, 10}}),
                                         anchor_of(2, "op2.4.SA", true, {{0, 1, 0, 9}})},
                                        1, writer);
        REQUIRE(rows == vector<string>{"op2.4.SA:1"});
        REQUIRE(writer.tally().anchors_left_out == 1);
    }
    SECTION("another alignment of the read's name, as a mate is, is another read") {
        const auto rows = write_anchors({anchor_of(3, "op3.0.SA", true, {{0, 2, 1, 10}}),
                                         anchor_of(2, "op2.4.SA", true, {{0, 1, 0, 9}})},
                                        1, writer);
        REQUIRE(rows == vector<string>{"op2.4.SA:1", "op3.0.SA:1"});
        REQUIRE(writer.tally().placements_left_out == 0);
    }
    SECTION("an off-panel anchor left below the minimum reads is left out, and frees its points") {
        const auto rows = write_anchors({anchor_of(2, "op2.4.SA", true, {{0, 1, 0, 9}, {1, 1, 0, 30}}),
                                         anchor_of(3, "op3.0.SA", true, {{1, 1, 0, 30}, {2, 1, 0, 40}}),
                                         anchor_of(1, ">1>2", false, {{0, 1, 0, 9}, {3, 1, 0, 60}})},
                                        2, writer);
        REQUIRE(rows == vector<string>{">1>2:2", "op3.0.SA:2"});
        REQUIRE(writer.tally().placements_left_out == 1);
        REQUIRE(writer.tally().anchors_left_out == 1);
    }
}

}
}
