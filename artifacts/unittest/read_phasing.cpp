/// \file unittest/read_phasing.cpp
///
/// Unit tests for read phasing, which decides from the reads which heterozygous sites to flip.
///
/// Everything here is parity arithmetic over a chain: a sign error in one link flips every site
/// after it, and nothing in the VCF looks unusual. These tests pin the convention.
///

#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../read_phasing.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A site whose reads all carry slot `which` -- perfectly discriminating, never mismapped, so the
/// link either side of it is as strong as the model allows.
static PhaseSite site(size_t key, size_t pos, int which, size_t n_reads = 20,
                      double reliability = 12.0, uint64_t first_read = 0) {
    PhaseSite s;
    s.record_key = key;
    s.phase_set = 1;
    s.position = pos;
    s.reliability = reliability;
    for (size_t i = 0; i < n_reads; ++i) {
        s.read_key.push_back(first_read + i);
        // Half the reads on each haplotype, so the site is a genuine het. `which` flips which
        // haplotype slot 0 names, which is the whole quantity under test.
        const bool upper = i >= n_reads / 2;
        const double q = upper ? 1.0 : 0.0;
        s.q0.push_back((float)(which ? 1.0 - q : q));
        s.c.push_back(0.95f);
    }
    return s;
}

TEST_CASE("read phasing leaves a chain the reads agree with alone", "[read_phasing]") {
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 0), site(3, 300, 0)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(flips.empty());
    REQUIRE(counters.reliable == 3);
    REQUIRE(counters.breaks == 0);
}

TEST_CASE("read phasing flips a site the reads put the other way round", "[read_phasing]") {
    // Site 2 is written with its slots swapped, so the reads say it is out of frame -- and so is
    // every site after it, because a switch propagates.
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 1), site(3, 300, 1)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(flips.count(2) == 1);
    REQUIRE(flips.count(3) == 1);
    REQUIRE(flips.count(1) == 0);
    REQUIRE(counters.breaks == 0);
}

TEST_CASE("the chain steps over an unreliable site rather than through it", "[read_phasing]") {
    // THE point of the design. The middle site is below the confidence threshold, so it must not
    // carry the link: the outer two are joined directly and the middle one is hung off the result.
    // Its own reads are still consistent here, so it should land the right way up -- but a mistake
    // on it is a flip, not a switch, and the outer relation survives either way.
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 0, 20, 1.0), site(3, 300, 0)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(counters.reliable == 2);
    REQUIRE(counters.hung == 1);
    REQUIRE(flips.empty());
}

TEST_CASE("a site sharing no reads with the chain cannot move it", "[read_phasing]") {
    // Disjoint read keys, so there is no evidence at all. The panel's frame has to stand: returning
    // a flip here would invent a phase out of nothing.
    vector<PhaseSite> sites{site(1, 100, 0, 20, 12.0, 0), site(2, 200, 1, 20, 12.0, 500)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(flips.empty());
    REQUIRE(counters.breaks == 1);
    REQUIRE(counters.breaks_no_reads == 1);
}

TEST_CASE("phase is not compared across blocks", "[read_phasing]") {
    // Two sites in different phase sets. Whatever the reads say, there is no relation to decide:
    // each block's first site pins its own frame.
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 1)};
    sites[1].phase_set = 2;
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(flips.empty());
    REQUIRE(counters.chains == 2);
}

TEST_CASE("break counts describe the chain the last coherence round left", "[read_phasing]") {
    // Two phase sets built alike. In each, A, B and D share reads that agree, C's reads alternate
    // between its alleles, so the coherence step removes it, and E shares no reads. The first pass
    // breaks the chain on both sides of C and before E; the second, without C, breaks only
    // before E, with no read across that break.
    vector<PhaseSite> sites;
    for (size_t ps : {1, 2}) {
        const uint64_t reads = ps * 1000;
        const size_t key = ps * 10;
        PhaseSite c = site(key + 3, 300, 0, 20, 12.0, reads);
        for (size_t i = 0; i < c.q0.size(); ++i) {
            c.q0[i] = (i % 2 == 0) ? 1.0f : 0.0f;
        }
        sites.push_back(site(key + 1, 100, 0, 20, 12.0, reads));
        sites.push_back(site(key + 2, 200, 0, 20, 12.0, reads));
        sites.push_back(c);
        sites.push_back(site(key + 4, 400, 0, 20, 12.0, reads));
        sites.push_back(site(key + 5, 500, 0, 20, 12.0, reads + 500));
        for (size_t i = sites.size() - 5; i < sites.size(); ++i) {
            sites[i].phase_set = ps;
        }
    }
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    REQUIRE(counters.chains == 2);
    REQUIRE(counters.demoted_incoherent == 2);
    REQUIRE(counters.breaks == 2);
    REQUIRE(counters.breaks_no_reads == 2);
    // One round in each chain.
    REQUIRE(counters.coherence_rounds_run == 1);
}

TEST_CASE("a hung site that no read links keeps the panel's order", "[read_phasing]") {
    // With `panel_weight` at 0 and no shared read, the hang stage has no vote either way, and the
    // panel gave the site an order. Flipping it would invent a phase out of nothing.
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 1, 20, 1.0, 500), site(3, 300, 0)};
    ReadPhasingParams params;
    params.panel_weight = 0.0;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(counters.hung == 1);
    REQUIRE(counters.hung_no_reads == 1);
    REQUIRE(flips.empty());
}

TEST_CASE("reads order the sites the panel cannot that the chain does not reach", "[read_phasing]") {
    // Chain sites 1 and 3 share no read with sites 10-13, which no panel haplotype orders and which
    // are too unreliable to join the chain. 13 shares reads with 4, which is hung from the chain by
    // the panel's vote, and takes its order from it. 10 and 11 share reads only with each other, so
    // they are ordered against each other, 11 the other way round from 10. 12 shares a read with
    // nothing, and stays unordered.
    vector<PhaseSite> sites{site(1, 100, 0), site(3, 300, 0),
                            site(4, 340, 0, 20, 1.0, 1000),
                            site(10, 150, 0, 20, 1.0, 500), site(11, 250, 1, 20, 1.0, 500),
                            site(12, 280, 0, 20, 1.0, 900), site(13, 320, 1, 20, 1.0, 1000)};
    for (PhaseSite& s : sites) {
        s.panel_ordered = s.record_key < 10;
    }
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    unordered_set<size_t> read_ordered;
    const auto flips = read_phase_flips(sites, params, counters, &read_ordered);
    REQUIRE(counters.reliable == 2);
    REQUIRE(counters.hung == 5);
    REQUIRE(counters.hung_no_reads == 2);
    REQUIRE(read_ordered == unordered_set<size_t>{10, 11, 13});
    REQUIRE(flips == unordered_set<size_t>{11, 13});
}

TEST_CASE("a run of the chain takes the panel's order from its first site the panel orders", "[read_phasing]") {
    // Two runs no read joins. Each starts with a reliable site the panel does not order, 10 and 11,
    // whose reads put it the other way round from the panel sites after it. The panel sites keep
    // the panel's order, and 10 and 11 are flipped to match them and count as ordered by reads.
    vector<PhaseSite> sites{site(10, 100, 1), site(2, 200, 0), site(3, 300, 0),
                            site(11, 400, 1, 20, 12.0, 500), site(4, 500, 0, 20, 12.0, 500)};
    for (PhaseSite& s : sites) {
        s.panel_ordered = s.record_key < 10;
    }
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    unordered_set<size_t> read_ordered;
    const auto flips = read_phase_flips(sites, params, counters, &read_ordered);
    REQUIRE(counters.reliable == 5);
    REQUIRE(counters.breaks_no_reads == 1);
    REQUIRE(read_ordered == unordered_set<size_t>{10, 11});
    REQUIRE(flips == unordered_set<size_t>{10, 11});
}

TEST_CASE("each site takes the run of the chain its order is relative to", "[read_phasing]") {
    // Sites 1-3 share reads 0-19 and sites 4-5 reads 500-519, so no read joins the two runs. Site
    // 6 is unreliable and shares reads only with the second run, though the nearest site the panel
    // orders is in the first, so the reads that order it vote from the second run. Site 7 shares no
    // read with anything and keeps the panel's order relative to site 1, the nearest site the
    // panel orders.
    vector sites{site(1, 100, 0), site(7, 120, 0, 20, 1.0, 900), site(2, 200, 0),
                            site(3, 300, 0), site(6, 350, 0, 20, 1.0, 500),
                            site(4, 400, 0, 20, 12.0, 500), site(5, 500, 0, 20, 12.0, 500)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    REQUIRE(counters.breaks_no_reads == 1);
    REQUIRE(counters.hung == 2);
    unordered_map<size_t, size_t> run;
    for (const PhaseSite& s : sites) {
        run[s.record_key] = s.run;
    }
    REQUIRE(run == unordered_map<size_t, size_t>{{1, 0}, {2, 0}, {3, 0}, {7, 0},
                                                 {4, 1}, {5, 1}, {6, 1}});
}

TEST_CASE("a hung site keeps the run of the reads that agree with it, however weak",
          "[read_phasing]") {
    // Site 6 has one read, which it shares with sites 4 and 5 of the second run and which agrees
    // with them. The panel's vote, toward site 3 of the first run, outweighs its read's vote from
    // the second run but agrees with it, so the site's order is the one its read gives it.
    vector sites{site(1, 100, 0), site(2, 200, 0), site(3, 300, 0),
                            site(6, 350, 1, 1, 1.0, 500), site(4, 400, 0, 20, 12.0, 500),
                            site(5, 500, 0, 20, 12.0, 500)};
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(flips.empty());
    for (const PhaseSite& s : sites) {
        REQUIRE(s.run == (s.record_key <= 3 ? 0 : 1));
    }
}

TEST_CASE("a read votes once on a hung site however many of one run's phase-chain sites it spans",
          "[read_phasing]") {
    // Six reliable sites in one run whose 20 reads all span them, and an unreliable site in their
    // middle sharing one of those reads, read 0, which puts it the other way round from the panel.
    // Its link to each of the six is about -1.3, so the six links would outvote a panel vote of 3.
    // As one read its vote is about -1.6, so the site keeps the panel's order against a vote of 3
    // and takes the read's against a vote of 1.
    vector<PhaseSite> sites;
    for (size_t k = 1; k <= 6; ++k) {
        sites.push_back(site(k, 100 * k, 0));
    }
    sites.insert(sites.begin() + 3, site(7, 350, 0, 1, 1.0));
    REQUIRE(phase_link(sites[3], sites[0], 0.0) < -1.0);
    ReadPhasingParams params;
    REQUIRE(6 * phase_link(sites[3], sites[0], 0.0) + params.panel_weight < 0.0);
    vector<PhaseSite> again = sites;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(!counters.panel.fitted);
    REQUIRE(counters.hung == 1);
    REQUIRE(flips.empty());
    params.panel_weight = 1.0;
    REQUIRE(read_phase_flips(again, params, counters) == unordered_set<size_t>{7});
}

TEST_CASE("a read that cannot tell a hung site's alleles apart does not link it to the chain",
          "[read_phasing]") {
    // Site 7, which no panel haplotype orders, shares read 0 with the chain, but the read fits its
    // two alleles equally well. The site stays unordered rather than taking an arbitrary order
    // from that read.
    vector sites{site(1, 100, 0), site(2, 200, 0), site(7, 250, 0, 1, 1.0), site(3, 300, 0)};
    sites[2].panel_ordered = false;
    sites[2].q0[0] = 0.5f;
    REQUIRE(phase_link(sites[2], sites[0], 0.0) == 0.0);
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    unordered_set<size_t> read_ordered;
    const auto flips = read_phase_flips(sites, params, counters, &read_ordered);
    REQUIRE(counters.hung == 1);
    REQUIRE(counters.hung_no_reads == 1);
    REQUIRE(!read_ordered.contains(7));
    REQUIRE(flips.empty());
}

TEST_CASE("sites only reads among themselves order make a run of their own", "[read_phasing]") {
    // Sites 8 and 9 are unreliable, the panel does not order them, and they share reads only with
    // each other, so they are ordered against each other and are a run apart from the chain's.
    vector sites{site(1, 100, 0), site(2, 200, 0), site(8, 300, 0, 20, 1.0, 700),
                            site(9, 400, 0, 20, 1.0, 700)};
    sites[2].panel_ordered = false;
    sites[3].panel_ordered = false;
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    unordered_map<size_t, size_t> run;
    for (const PhaseSite& s : sites) {
        run[s.record_key] = s.run;
    }
    REQUIRE(run == unordered_map<size_t, size_t>{{1, 0}, {2, 0}, {8, 1}, {9, 1}});
}

/// `count` sites 1 kb apart from position `start`, each sharing 20 reads with the next, so that
/// reads join them into one piece. The panel orders site k's slots the wrong way round where
/// `swapped(k)`. Even read keys, counted from `first_read`, come from strand 0 and odd ones from
/// strand 1.
template <typename Swapped>
static vector<PhaseSite> linked_sites(size_t count, size_t start, uint64_t first_read,
                                      const Swapped& swapped) {
    vector<PhaseSite> sites;
    for (size_t k = 0; k < count; ++k) {
        PhaseSite s;
        s.record_key = start / 1000 + k + 1;
        s.phase_set = 1;
        s.position = start + 1000 * k;
        s.reliability = 12.0;
        for (uint64_t read = first_read + 20 * k; read < first_read + 20 * k + 40; ++read) {
            s.read_key.push_back(read);
            s.q0.push_back((read % 2 == 1) == swapped(k) ? 1.0f : 0.0f);
            s.c.push_back(0.95f);
        }
        sites.push_back(s);
    }
    return sites;
}

TEST_CASE("read phasing fits the panel's error rates from its pieces", "[read_phasing]") {
    // 399 neighbouring pairs 1 kb apart, the panel switching between 9 of them.
    vector<PhaseSite> sites =
        linked_sites(400, 1000, 0, [](size_t k) { return (k / 40) % 2 == 1; });
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    REQUIRE(counters.breaks == 0);
    REQUIRE(counters.panel_pairs == 399);
    REQUIRE(counters.panel_pairs_reversed == 9);
    REQUIRE(counters.panel.fitted);
    // Every pair is 1 kb apart, so the fit sees only the agreement at 1 kb.
    const double agreement = (1.0 + counters.panel.correlation(1000)) / 2.0;
    REQUIRE(agreement == Approx(390.0 / 399.0).epsilon(0.003));
    REQUIRE(counters.panel_fit.size() == 8);

    // Too few pairs fit no rates.
    vector<PhaseSite> few = linked_sites(100, 1000, 0, [](size_t) { return false; });
    ReadPhasingCounters few_counters;
    read_phase_flips(few, params, few_counters);
    REQUIRE(few_counters.panel_pairs == 99);
    REQUIRE(!few_counters.panel.fitted);

    // Pairs the panel never reverses still leave it a flip rate, so that its vote stays finite.
    vector<PhaseSite> clean = linked_sites(300, 1000, 0, [](size_t) { return false; });
    ReadPhasingCounters clean_counters;
    read_phase_flips(clean, params, clean_counters);
    REQUIRE(clean_counters.panel.fitted);
    REQUIRE(clean_counters.panel.flip_rate > 0.0);
    REQUIRE(std::isfinite(clean_counters.panel.vote(0)));
}

TEST_CASE("the panel orients a run against the one before it at the break between them",
          "[read_phasing]") {
    // Two runs no read joins, the second 1 kb after the first. The panel switches halfway through
    // the first run and has no error from there on. Each run's own first site would say the
    // second run keeps the panel's order; the sites either side of the break say that it is flipped
    // with the end of the first run.
    auto swapped = [](size_t k) { return k >= 150; };
    vector<PhaseSite> sites = linked_sites(300, 1000, 0, swapped);
    vector<PhaseSite> second =
        linked_sites(300, 301000, 1000000, [](size_t) { return true; });
    sites.insert(sites.end(), second.begin(), second.end());
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(counters.panel.fitted);
    REQUIRE(counters.breaks_no_reads == 1);
    for (const PhaseSite& s : sites) {
        // Sites from the 151st on are the other way round from the panel, in both runs.
        REQUIRE(flips.contains(s.record_key) == (s.position >= 151000));
        REQUIRE(s.run == (s.position > 300000 ? 1 : 0));
        if (s.run == 1) {
            REQUIRE(s.run_correlation == Approx(counters.panel.correlation(1000)));
        }
    }
}


TEST_CASE("a break the panel decides against the reads starts a run", "[read_phasing]") {
    // Two pieces that one weak read joins, saying their sites either side of the break are the
    // other way round from the panel's order. The panel's vote outweighs it, so the later piece
    // keeps the panel's order, and starts a run whose correlation says how likely that is right.
    vector<PhaseSite> sites = linked_sites(300, 1000, 0, [](size_t) { return false; });
    vector<PhaseSite> second = linked_sites(300, 301000, 1000000, [](size_t) { return false; });
    sites.insert(sites.end(), second.begin(), second.end());
    const uint64_t weak = 5000000;
    for (size_t k : {size_t(299), size_t(300)}) {
        sites[k].read_key.push_back(weak);
        sites[k].q0.push_back(k == 299 ? 1.0f : 0.0f);
        sites[k].c.push_back(0.6f);
    }
    REQUIRE(phase_link(sites[299], sites[300], 0.0) < 0.0);
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    const auto flips = read_phase_flips(sites, params, counters);
    REQUIRE(counters.panel.fitted);
    REQUIRE(counters.breaks == 1);
    REQUIRE(counters.breaks_no_reads == 0);
    REQUIRE(flips.empty());
    for (const PhaseSite& s : sites) {
        REQUIRE(s.run == (s.position > 300000 ? 1 : 0));
        if (s.run == 1) {
            REQUIRE(s.run_correlation > 0.99);
            REQUIRE(s.run_correlation < 1.0);
        }
    }
}

TEST_CASE("the panel relates a run only to sites oriented against the run before it",
          "[read_phasing]") {
    // Three runs no read joins. The panel orders no site of the middle one, so it relates neither
    // neighbour to the middle run, and its order of the sites either side of the middle run says
    // nothing about how either is oriented against it.
    vector<PhaseSite> sites = linked_sites(300, 1000, 0, [](size_t) { return false; });
    vector<PhaseSite> middle = linked_sites(2, 301000, 1000000, [](size_t) { return false; });
    for (PhaseSite& s : middle) {
        s.panel_ordered = false;
    }
    vector<PhaseSite> last = linked_sites(300, 303000, 2000000, [](size_t) { return true; });
    sites.insert(sites.end(), middle.begin(), middle.end());
    sites.insert(sites.end(), last.begin(), last.end());
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    REQUIRE(counters.panel.fitted);
    REQUIRE(counters.breaks_no_reads == 2);
    for (const PhaseSite& s : sites) {
        // The run starts at 1000, 301000 and 303000.
        const size_t run = (s.position >= 301000 ? 1 : 0) + (s.position >= 303000 ? 1 : 0);
        REQUIRE(s.run == run);
        REQUIRE(s.run_correlation == 0.0);
    }
}

TEST_CASE("a chain site the panel does not order is ordered by reads only with company", "[read_phasing]") {
    // 20 and 21, which the panel does not order, form a run of their own, so reads order 21
    // against 20, but nothing orders the run as a whole. 22, alone in a run, is not ordered at all.
    // 5, hung with no reads, takes the panel's vote relative to 2, the nearest chain site the
    // panel orders, not relative to 21, which is nearer but whose order the panel did not give.
    vector<PhaseSite> sites{site(2, 100, 0), site(20, 300, 0, 20, 12.0, 900),
                            site(21, 310, 1, 20, 12.0, 900), site(5, 320, 0, 20, 1.0, 2000),
                            site(22, 400, 0, 20, 12.0, 3000)};
    for (PhaseSite& s : sites) {
        s.panel_ordered = s.record_key < 10;
    }
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    unordered_set<size_t> read_ordered;
    const auto flips = read_phase_flips(sites, params, counters, &read_ordered);
    REQUIRE(counters.hung == 1);
    REQUIRE(counters.hung_no_reads == 1);
    REQUIRE(read_ordered == unordered_set<size_t>{20, 21});
    REQUIRE(flips == unordered_set<size_t>{21});
}

TEST_CASE("paired mates at one site are one read, the more confident mate", "[read_phasing]") {
    // Mates share a read key. Read 5 has a confident mate and a weak one that disagrees with it.
    PhaseSite forward;
    forward.read_key = {5, 6, 5};
    forward.q0 = {0.9f, 0.5f, 0.2f};
    forward.c = {0.95f, 0.9f, 0.5f};
    PhaseSite backward;
    backward.read_key = {5, 6, 5};
    backward.q0 = {0.2f, 0.5f, 0.9f};
    backward.c = {0.5f, 0.9f, 0.95f};
    for (PhaseSite* s : {&forward, &backward}) {
        merge_mates(*s);
        REQUIRE(s->read_key == vector<uint64_t>{5, 6});
        REQUIRE(s->q0 == vector<float>{0.9f, 0.5f});
        REQUIRE(s->c == vector<float>{0.95f, 0.9f});
    }

    // The site's reliability averages reads, not mates: two rows of read 1, one clean and one that
    // fits both alleles equally, count as the clean one alone.
    PhaseReadEvidence ev;
    ev.n_alleles = 2;
    ev.length_weighted = false;
    ev.read_key = {1, 1, 2};
    ev.mismap = {0.02f, 0.02f, 0.02f};
    ev.rel = {1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f};
    const PhaseSite reduced = reduce_to_pair(ev, 0, 1);
    REQUIRE(reduced.read_key == vector<uint64_t>{1, 2});
    REQUIRE(reduced.q0[0] == 1.0f);
    REQUIRE(reduced.reliability == Approx(read_confidence(reduced.q0[0], reduced.c[0])));
    REQUIRE(read_confidence(reduced.q0[0], reduced.c[0])
            == Approx(read_confidence(reduced.q0[1], reduced.c[1])));

    // Read phasing merges the mates of a site built by hand too, so its links and coherence see
    // each read once.
    vector<PhaseSite> sites{site(1, 100, 0), site(2, 200, 0), site(3, 300, 0)};
    sites[1].read_key.push_back(3);
    sites[1].q0.push_back(0.5f);
    sites[1].c.push_back(0.5f);
    ReadPhasingParams params;
    ReadPhasingCounters counters;
    read_phase_flips(sites, params, counters);
    REQUIRE(sites[1].read_key.size() == 20);
    REQUIRE(sites[1].q0[3] == 0.0f);
}

TEST_CASE("a link's sign is what the log odds say", "[read_phasing]") {
    const PhaseSite a = site(1, 100, 0);
    const PhaseSite same = site(2, 200, 0);
    const PhaseSite other = site(3, 300, 1);
    REQUIRE(phase_link(a, same, 0.0) > 0.0);
    REQUIRE(phase_link(a, other, 0.0) < 0.0);
    // Symmetric, and capped where asked.
    REQUIRE(phase_link(a, same, 0.0) == Approx(phase_link(same, a, 0.0)));
    REQUIRE(phase_link(a, same, 1.0) == Approx(1.0));
    REQUIRE(phase_link(a, other, 1.0) == Approx(-1.0));
}

}
}
