/// \file unittest/regenotype.cpp
///
/// Unit tests for the phase-aware genotype correction.
///
/// A run with a temper of 0 does not exercise `Lambda`, since at tau = 0 the weighting is the
/// identity whatever `Lambda` holds. A run limited to one linkage pass does not either, since it
/// computes and reports the correction but does not apply it.
///
/// At a heterozygous genotype, the correction uses whichever order of the two alleles on the
/// strands fits the reads better. So a `Lambda` whose sign is wrong at every site gives the same
/// corrected likelihoods as the right one. A `Lambda` wired backwards at some sites would look
/// like a run with `RegenotypeParams::shuffle` on. These tests are what check the sign of
/// `Lambda`.

#include <cmath>
#include <map>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../anchor.hpp"
#include "../read_strand_table.hpp"
#include "../regenotype.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A site's retained evidence: `n_reads` reads over two alleles, read `i` carrying allele
/// `i % 2` cleanly. `len0`/`len1` are the alleles' spelled lengths, which is what decides whether
/// the slot weights are equal.
static PhaseReadEvidence evidence(size_t n_reads, uint32_t len0, uint32_t len1,
                                  float mismap = 0.02f, float mean_read_length = 15000.0f) {
    PhaseReadEvidence ev;
    ev.n_alleles = 2;
    ev.allele_length = {len0, len1};
    ev.mean_read_length = mean_read_length;
    ev.length_weighted = true;
    for (size_t i = 0; i < n_reads; ++i) {
        ev.read_key.push_back((uint64_t)(i + 1));
        ev.mismap.push_back(mismap);
        const bool carries_1 = (i % 2) == 1;
        ev.rel.push_back(carries_1 ? 0.0f : 1.0f);
        ev.rel.push_back(carries_1 ? 1.0f : 0.0f);
    }
    return ev;
}

/// Every read confidently on the strand matching the allele it carries, from other sites.
static LambdaTable confident(const PhaseReadEvidence& ev, double magnitude = 4.0, size_t sites = 8) {
    LambdaTable t;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = (i % 2 == 0) ? magnitude : -magnitude;
        rl.sites = sites;
        rl.phase_set = 1;
        t[ev.read_key[i]] = rl;
    }
    return t;
}

static map<vector<int>, double> flat_gl() {
    return {{{0, 0}, -50.0}, {{0, 1}, -20.0}, {{1, 1}, -50.0}};
}

TEST_CASE("a temper of zero leaves every genotype likelihood bit-identical",
          "[regenotype]") {
    // At tau = 0 the tilted weights equal the slot weights, so the corrected and uncorrected
    // mixtures are the same arithmetic on the same values. The likelihoods must therefore be
    // equal, not merely close.
    RegenotypeParams params;
    RegenotypeCounters counters;
    for (uint32_t len1 : {10u, 400u, 12000u}) {
        PhaseReadEvidence ev = evidence(30, 10, len1);
        LambdaTable lambda = confident(ev);
        map<vector<int>, double> gl = flat_gl();
        const map<vector<int>, double> before = gl;
        phase_aware_correction(ev, lambda, 1, NO_RUN, 0, {},
                               StrandCalibration{0.0, 1.0}, params, gl, counters);
        for (const auto& kv : before) {
            REQUIRE(gl.at(kv.first) == kv.second);
        }
    }
    // The two orders of the alleles on the strands must also score equal, so that no site
    // counts as preferring the reversed order. A difference in the last bit would make the
    // correction at temper 0 near zero rather than zero, and a temper of 0 must change nothing.
    REQUIRE(counters.order_reversed == 0);
}

TEST_CASE("homozygous likelihoods never move, at any temper", "[regenotype]") {
    // A homozygote's mixture collapses to rel(r, a) whatever the weights are. This is what
    // confines the correction to re-ranking hets against each other and hets against homs.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 12000);
    LambdaTable lambda = confident(ev);
    for (double tau : {0.25, 1.0, 50.0}) {
        map<vector<int>, double> gl = flat_gl();
        phase_aware_correction(ev, lambda, 1, NO_RUN, 0, {},
                               StrandCalibration{tau, 1.0}, params, gl, counters);
        REQUIRE(gl.at({0, 0}) == -50.0);
        REQUIRE(gl.at({1, 1}) == -50.0);
        REQUIRE(gl.at({0, 1}) > -20.0);   // and the het does move
    }
}

TEST_CASE("a read that spans nothing else contributes exactly nothing", "[regenotype]") {
    // Leave-one-out gives such a read Lambda = 0, at which point its tilted weights ARE the slot
    // weights and its term cancels between the two sides of the delta. This is what makes the
    // correction degenerate harmlessly on short reads, where most reads span one site.
    //
    // The length ratio is deliberately extreme: it is the case where pinning the weights to the
    // SLOT rather than the ALLELE would break the property, because the two orders would then
    // score differently on reads carrying no phase information at all.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 12000);
    LambdaTable lambda;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = 0.0;
        rl.sites = 1;
        rl.phase_set = 1;
        lambda[ev.read_key[i]] = rl;
    }
    map<vector<int>, double> gl = flat_gl();
    const map<vector<int>, double> before = gl;
    phase_aware_correction(ev, lambda, 1, NO_RUN, 0, {}, StrandCalibration{2.0, 1.0}, params, gl,
                           counters);
    for (const auto& kv : before) {
        REQUIRE(gl.at(kv.first) == kv.second);
    }
}

TEST_CASE("reversing the pair leaves the uncorrected mixture bit-identical", "[regenotype]") {
    // The property above rests on `allele_length_weights` deriving each slot's weight from that slot's
    // allele, so reordering the pair reorders the weights with it. Keeping `w = (w0, w1)` across the
    // swap would make a read that spans no other site depend on the order, with no phase
    // information in it, even at tau = 0.
    PhaseReadEvidence ev = evidence(6, 10, 12000);
    const vector<double> fwd = allele_length_weights(ev.allele_length, ev.n_alleles,
                                                 ev.mean_read_length, ev.length_weighted,
                                                 vector<int>{0, 1});
    const vector<double> rev = allele_length_weights(ev.allele_length, ev.n_alleles,
                                                 ev.mean_read_length, ev.length_weighted,
                                                 vector<int>{1, 0});
    REQUIRE(fwd.size() == 2);
    REQUIRE(rev.size() == 2);
    REQUIRE(fwd[0] == rev[1]);
    REQUIRE(fwd[1] == rev[0]);
    // And they are genuinely unequal here, so the test can fail.
    REQUIRE(fwd[0] != fwd[1]);
}

TEST_CASE("leave-one-out equals accumulating without the site", "[regenotype]") {
    // The subtraction is the whole of the leave-one-out, and it has to be exact rather than
    // close: if a site's own evidence leaks into the prior that re-genotypes it, the site
    // confirms itself and the arm measures its own assumption.
    vector<PhaseSite> sites;
    for (size_t s = 0; s < 4; ++s) {
        PhaseSite site;
        site.record_key = s;
        site.phase_set = 1;
        site.position = 100 * s;
        for (size_t i = 0; i < 10; ++i) {
            site.read_key.push_back((uint64_t)(i + 1));
            site.q0.push_back((i % 2 == 0) ? 0.98f : 0.02f);
            site.c.push_back(0.9f);
        }
        sites.push_back(site);
    }
    RegenotypeCounters counters;
    LambdaTable all;
    accumulate_lambda(sites, {}, all, counters);

    // Drop site 2 and accumulate again: that is what the subtraction must reproduce.
    vector<PhaseSite> without;
    for (const PhaseSite& s : sites) {
        if (s.record_key != 2) {
            without.push_back(s);
        }
    }
    RegenotypeCounters c2;
    LambdaTable minus_one;
    accumulate_lambda(without, {}, minus_one, c2);

    unordered_map<uint64_t, double> own;
    site_own_log_odds(sites[2], false, own);
    for (const auto& kv : minus_one) {
        const double by_subtraction = all.at(kv.first).lambda - own.at(kv.first);
        REQUIRE(by_subtraction == Approx(kv.second.lambda).epsilon(1e-12));
    }
}

TEST_CASE("a flipped site enters Lambda with the opposite sign", "[regenotype]") {
    // `read_phase_flips` swaps a site's slot order after the PhaseSites were built, so a site in the
    // flip set describes the panel's order, not the chosen one, and its contribution must enter
    // with the opposite sign.
    PhaseSite site;
    site.record_key = 7;
    site.phase_set = 1;
    site.position = 0;
    for (size_t i = 0; i < 8; ++i) {
        site.read_key.push_back((uint64_t)(i + 1));
        site.q0.push_back(0.95f);
        site.c.push_back(0.9f);
    }
    RegenotypeCounters ca, cb;
    LambdaTable plain, flipped;
    accumulate_lambda({site}, {}, plain, ca);
    accumulate_lambda({site}, {site.record_key}, flipped, cb);
    for (const auto& kv : plain) {
        REQUIRE(flipped.at(kv.first).lambda == Approx(-kv.second.lambda).epsilon(1e-12));
    }
    REQUIRE(plain.at(1).lambda > 0.0);
}

TEST_CASE("a paired mate counts once", "[regenotype]") {
    // Paired mates share a read name and so a read_key, but the read source keeps both, as separate
    // rows under one key. A fragment lies on one strand, so it counts once. `merge_mates` leaves
    // one row per key, and a site that still holds several adds the key once.
    PhaseSite site;
    site.record_key = 1;
    site.phase_set = 1;
    site.position = 0;
    for (size_t i = 0; i < 4; ++i) {
        site.read_key.push_back(42);   // one fragment, four rows
        site.q0.push_back(0.95f);
        site.c.push_back(0.9f);
    }
    PhaseSite single = site;
    single.read_key.assign(1, 42);
    single.q0.assign(1, 0.95f);
    single.c.assign(1, 0.9f);

    RegenotypeCounters ca, cb;
    LambdaTable many, one;
    accumulate_lambda({site}, {}, many, ca);
    accumulate_lambda({single}, {}, one, cb);
    REQUIRE(many.size() == 1);
    REQUIRE(many.at(42).lambda == Approx(one.at(42).lambda).epsilon(1e-12));
    REQUIRE(many.at(42).sites == 1);
}

TEST_CASE("the correction rewards a phase-coherent split and punishes an incoherent one",
          "[regenotype]") {
    // The direction of the whole thing, and it is not the obvious one. A het whose reads split
    // along the phase gains; a het whose reads do not is pushed toward hom. So the exposure is
    // recall rather than precision, and the het count should be expected to fall.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 10);   // equal lengths, so the weights are 1/2 each

    LambdaTable coherent = confident(ev);
    map<vector<int>, double> gl_coherent = flat_gl();
    phase_aware_correction(ev, coherent, 1, NO_RUN, 0, {},
                           StrandCalibration{1.0, 1.0}, params, gl_coherent, counters);

    // Same magnitudes, strand assignment uncorrelated with the allele carried.
    LambdaTable incoherent = confident(ev);
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        incoherent[ev.read_key[i]].lambda = ((i / 2) % 2 == 0) ? 4.0 : -4.0;
    }
    map<vector<int>, double> gl_incoherent = flat_gl();
    phase_aware_correction(ev, incoherent, 1, NO_RUN, 0, {},
                           StrandCalibration{1.0, 1.0}, params, gl_incoherent, counters);

    const double het_gap_coherent = gl_coherent.at({0, 1}) - gl_coherent.at({0, 0});
    const double het_gap_incoherent = gl_incoherent.at({0, 1}) - gl_incoherent.at({0, 0});
    const double het_gap_before = -20.0 - -50.0;
    REQUIRE(het_gap_coherent > het_gap_before);
    REQUIRE(het_gap_incoherent < het_gap_coherent);
}

TEST_CASE("at temper 0 every read's log-odds are 0, whatever the ceiling", "[regenotype]") {
    // A temper of 0 must change nothing, so every read's log-odds must be 0 there for any
    // ceiling c. At tau = 0 the sigmoid is 1/2, so the probability is c / 2 + (1 - c) / 2 = 1/2,
    // and its logit is 0.
    for (double c : {1.0, 0.99, 0.95, 0.5, 0.01}) {
        for (double lam : {0.0, 17.3, -1215.0, 1e9, -1e-9}) {
            REQUIRE(StrandCalibration{0.0, c}.log_odds(lam) == 0.0);
        }
    }
}

TEST_CASE("the ceiling caps how confident a read may be", "[regenotype]") {
    // A read with |Lambda| = 1215, at a temper of about 0.05, has tau * Lambda of about 60. With
    // a ceiling of 1, its strand probability is within 1e-12 of 1, so its log-odds are large.
    // With a ceiling of 0.95, the probability is at most 0.975, and the log-odds are at most
    // about 3.7.
    const double tau = 0.0498438;
    const double wild = StrandCalibration{tau, 1.0}.log_odds(1215.0);
    const double held = StrandCalibration{tau, 0.95}.log_odds(1215.0);
    REQUIRE(wild > 20.0);
    REQUIRE(held < 5.0);
    REQUIRE(std::isfinite(wild));   // clamped off 1, or the logit is infinite
    // Under the ceiling, the log-odds still grow with |Lambda|, and negating Lambda negates them.
    REQUIRE(StrandCalibration{tau, 0.95}.log_odds(102.0) < held);
    REQUIRE(StrandCalibration{tau, 0.95}.log_odds(-1215.0) == Approx(-held).epsilon(1e-12));
}

TEST_CASE("the haploid inclusion weight is the identity at temper 0", "[regenotype]") {
    // The parameterisation exists for this. Written as the read's POSTERIOR for this strand the
    // inclusion would be 1/2 at tau = 0 and half of every read's weight would vanish on a run
    // that is meant to change nothing; capped at 1 it is exactly 1 there.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 12000);
    LambdaTable lambda = confident(ev);
    map<vector<int>, double> gl = {{{0}, -20.0}, {{1}, -35.0}};
    const map<vector<int>, double> before = gl;
    for (int sign : {1, -1}) {
        gl = before;
        haploid_inclusion_correction(ev, lambda, 1, NO_RUN, {},
                                     StrandCalibration{0.0, 0.95}, sign, params, gl, counters);
        for (const auto& kv : before) {
            REQUIRE(gl.at(kv.first) == kv.second);
        }
    }
}

TEST_CASE("the haploid inclusion weight only removes evidence, never adds it", "[regenotype]") {
    // A read that fits the allele perfectly contributes the same whether it is included or not,
    // so excluding reads cannot manufacture support; a read that fits it not at all stops
    // penalising the allele as it is excluded, which is the entire mechanism.
    RegenotypeParams params;
    RegenotypeCounters counters;
    // ASYMMETRIC support -- 25 reads for allele 0, 5 for allele 1 -- because the property under
    // test is that excluding reads removes their DISCRIMINATING power. With the symmetric fixture
    // used elsewhere both alleles shift by the same amount and the gap is preserved, which says
    // nothing.
    PhaseReadEvidence ev;
    ev.n_alleles = 2;
    ev.allele_length = {10, 10};
    ev.mean_read_length = 15000.0f;
    ev.length_weighted = true;
    for (size_t i = 0; i < 30; ++i) {
        ev.read_key.push_back((uint64_t)(i + 1));
        ev.mismap.push_back(0.02f);
        const bool carries_1 = i >= 25;
        ev.rel.push_back(carries_1 ? 0.0f : 1.0f);
        ev.rel.push_back(carries_1 ? 1.0f : 0.0f);
    }
    // Every read placed on the OTHER strand, so every one of them is discounted.
    LambdaTable lambda;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = -60.0;
        rl.sites = 8;
        rl.phase_set = 1;
        lambda[ev.read_key[i]] = rl;
    }
    // The likelihoods must be the read term, as the direct pass's are: the correction subtracts the read
    // term, so likelihoods unrelated to `rel` would leave an arbitrary residue.
    auto read_term = [&](int a) {
        double t = 0.0;
        for (size_t r = 0; r < ev.num_reads(); ++r) {
            const double e = (double)ev.mismap[r];
            t += std::log((1.0 - e) * (double)ev.rel_at(r, (size_t)a) + e);
        }
        return t;
    };
    map<vector<int>, double> gl = {{{0}, read_term(0)}, {{1}, read_term(1)}};
    const double gap_before = gl.at({0}) - gl.at({1});
    REQUIRE(std::abs(gap_before) > 10.0);   // and the alleles really are discriminated to start
    haploid_inclusion_correction(ev, lambda, 1, NO_RUN, {},
                                 StrandCalibration{1.0, 1.0}, +1, params, gl, counters);
    // Both alleles rise, since reads that do not belong here no longer count against them. The
    // correction is added to the direct pass's likelihood, so the result is not bounded above by zero.
    REQUIRE(gl.at({0}) > read_term(0));
    REQUIRE(gl.at({1}) > read_term(1));
    // With every read excluded, every read term is log(1) = 0, and the two alleles are not told
    // apart at all.
    REQUIRE(std::abs(gl.at({0}) - gl.at({1})) < 1e-9);

    // Whereas with every read placed ON this strand, nothing is excluded and nothing moves.
    LambdaTable here;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = +60.0;
        rl.sites = 8;
        rl.phase_set = 1;
        here[ev.read_key[i]] = rl;
    }
    map<vector<int>, double> kept = {{{0}, read_term(0)}, {{1}, read_term(1)}};
    const map<vector<int>, double> kept_before = kept;
    haploid_inclusion_correction(ev, here, 1, NO_RUN, {},
                                 StrandCalibration{1.0, 1.0}, +1, params, kept, counters);
    REQUIRE(kept.at({0}) == kept_before.at({0}));
    REQUIRE(kept.at({1}) == kept_before.at({1}));
}

TEST_CASE("a read found in another phase set is not used", "[regenotype]") {
    // Each phase set labels its strands independently, so a read's strand from phase set 2 says
    // nothing about which strand of phase set 1 it is on.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable lambda;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = -60.0;
        rl.sites = 8;
        rl.phase_set = 2;
        lambda[ev.read_key[i]] = rl;
    }
    const map<vector<int>, double> before = {{{0}, -20.0}, {{1}, -35.0}};
    map<vector<int>, double> elsewhere = before;
    haploid_inclusion_correction(ev, lambda, 1, NO_RUN, {},
                                 StrandCalibration{1.0, 1.0}, +1, params, elsewhere, counters);
    REQUIRE(elsewhere == before);
    map<vector<int>, double> here = before;
    haploid_inclusion_correction(ev, lambda, 2, NO_RUN, {},
                                 StrandCalibration{1.0, 1.0}, +1, params, here, counters);
    REQUIRE(here != before);
    // A site with no phase set takes any read found in one phase set.
    map<vector<int>, double> unphased = before;
    haploid_inclusion_correction(ev, lambda, NO_PHASE_SET, NO_RUN, {},
                                 StrandCalibration{1.0, 1.0}, +1, params, unphased,
                                 counters);
    REQUIRE(unphased == here);

    map<vector<int>, double> gl = flat_gl();
    phase_aware_correction(ev, confident(ev), 2, NO_RUN, 0, {},
                           StrandCalibration{1.0, 1.0}, params, gl, counters);
    REQUIRE(gl == flat_gl());
}

TEST_CASE("a globally sign-flipped Lambda is invisible, and that is why these tests exist",
          "[regenotype]") {
    // Pinned deliberately. The max over the two orders absorbs a global sign flip, so no arm
    // downstream can detect one -- which is exactly why the leave-one-out and the flip-sign tests
    // above have to carry the weight.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable normal = confident(ev);
    LambdaTable flipped = confident(ev);
    for (auto& kv : flipped) {
        kv.second.lambda = -kv.second.lambda;
    }
    map<vector<int>, double> a = flat_gl(), b = flat_gl();
    phase_aware_correction(ev, normal, 1, NO_RUN, 0, {}, StrandCalibration{1.0, 1.0}, params, a,
                           counters);
    phase_aware_correction(ev, flipped, 1, NO_RUN, 0, {}, StrandCalibration{1.0, 1.0}, params, b,
                           counters);
    REQUIRE(a.at({0, 1}) == Approx(b.at({0, 1})).epsilon(1e-12));
}

TEST_CASE("the reads' order is compared with the chain's order, not the sorted one",
          "[regenotype]") {
    // In `confident`, reads carrying allele 0 are on strand 0, so the reads put allele 0 on
    // strand 0. They agree with a chain that does the same, and disagree with one that puts
    // allele 1 there. Negating Lambda moves allele 1 to strand 0, which then agrees with the
    // second chain.
    RegenotypeParams params;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable allele0_first = confident(ev);
    LambdaTable allele1_first = confident(ev);
    for (auto& kv : allele1_first) {
        kv.second.lambda = -kv.second.lambda;
    }
    auto reversed = [&](const LambdaTable& lambda, int strand0_allele) {
        RegenotypeCounters counters;
        map<vector<int>, double> gl = flat_gl();
        phase_aware_correction(ev, lambda, 1, NO_RUN, strand0_allele, {},
                               StrandCalibration{1.0, 1.0}, params, gl, counters);
        REQUIRE(counters.sites_corrected == 1);
        return counters.order_reversed;
    };
    REQUIRE(reversed(allele0_first, 0) == 0);
    REQUIRE(reversed(allele0_first, 1) == 1);
    REQUIRE(reversed(allele1_first, 1) == 0);
    REQUIRE(reversed(allele1_first, 0) == 1);
    // With no chain order there is nothing to disagree with.
    REQUIRE(reversed(allele1_first, -1) == 0);
}

TEST_CASE("the later rounds temper with the ceiling the first round fitted against",
          "[regenotype]") {
    // The first round fits the temper against the ceiling it was given, and the later rounds
    // reuse that fit rather than fitting again, so they must reuse the ceiling too.
    RegenotypeCounters first_round;
    first_round.calibration = StrandCalibration{0.05, 0.9};
    multipass::TemperFit fit;
    fit.keep(first_round);

    RegenotypeCounters later_round;
    fit.restore(later_round);
    REQUIRE(later_round.calibration.temper == 0.05);
    REQUIRE(later_round.calibration.ceiling == 0.9);
}

/// Put every read of `t` in run `run` of its phase set, all its sites there.
static LambdaTable in_run(LambdaTable t, size_t run) {
    for (auto& [key, rl] : t) {
        rl.runs[0] = run;
        rl.run_sites[0] = rl.sites;
        rl.run_lambda[0] = rl.lambda;
    }
    return t;
}

TEST_CASE("a site whose reads all lie in its own run is scored as if there were no runs",
          "[regenotype]") {
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 25);
    map<vector<int>, double> plain = flat_gl();
    phase_aware_correction(ev, confident(ev), 1, NO_RUN, 0, {}, StrandCalibration{1.0, 1.0},
                           params, plain, counters);
    map<vector<int>, double> own_run = flat_gl();
    phase_aware_correction(ev, in_run(confident(ev), 3), 1, 3, 0, {},
                           StrandCalibration{1.0, 1.0}, params, own_run, counters);
    // A site with no run of its own whose reads all lie in one run is one group too.
    map<vector<int>, double> no_site_run = flat_gl();
    phase_aware_correction(ev, in_run(confident(ev), 3), 1, NO_RUN, 0, {},
                           StrandCalibration{1.0, 1.0}, params, no_site_run, counters);
    REQUIRE(own_run == plain);
    REQUIRE(no_site_run == plain);
    REQUIRE(counters.sites_split_runs == 0);
}

TEST_CASE("reads from two runs are oriented apart", "[regenotype]") {
    // The second half of the reads lie in a run that the panel framed the other way round, so
    // their strands are all reversed. Scored as one group, every heterozygous order disagrees
    // with one half and the het loses. Scored per run, the first half takes the better order and
    // the second half's order is averaged over. Its reads all say the same, so it costs ln 2
    // against the halves' being in one frame.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable one_frame = confident(ev);
    LambdaTable inverted = confident(ev);
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda& rl = inverted[ev.read_key[i]];
        if (i >= 15) {
            rl.lambda = -rl.lambda;
        }
        rl.runs[0] = i >= 15 ? 1 : 0;
        rl.run_sites[0] = rl.sites;
        rl.run_lambda[0] = rl.lambda;
    }
    LambdaTable unlabelled = inverted;
    for (auto& [key, rl] : unlabelled) {
        rl.runs[0] = NO_RUN;
        rl.run_sites[0] = 0;
        rl.run_lambda[0] = 0.0;
    }
    const StrandCalibration calibration{1.0, 1.0};
    map<vector<int>, double> coherent = flat_gl();
    phase_aware_correction(ev, one_frame, 1, NO_RUN, 0, {}, calibration, params, coherent,
                           counters);
    map<vector<int>, double> one_group = flat_gl();
    phase_aware_correction(ev, unlabelled, 1, NO_RUN, 0, {}, calibration, params, one_group,
                           counters);
    REQUIRE(counters.sites_split_runs == 0);
    map<vector<int>, double> per_run = flat_gl();
    phase_aware_correction(ev, inverted, 1, 0, 0, {}, calibration, params, per_run, counters);
    REQUIRE(counters.sites_split_runs == 1);
    REQUIRE(one_group.at({0, 1}) < per_run.at({0, 1}));
    REQUIRE(per_run.at({0, 1}) == Approx(coherent.at({0, 1}) - log(2.0)));
    REQUIRE(per_run.at({0, 0}) == coherent.at({0, 0}));
    // At temper 0 the groups' sums cancel exactly.
    map<vector<int>, double> untempered = flat_gl();
    phase_aware_correction(ev, inverted, 1, 0, 0, {}, StrandCalibration{0.0, 1.0}, params,
                           untempered, counters);
    REQUIRE(untempered == flat_gl());
}

TEST_CASE("a read whose other sites lie in two runs keeps the larger part", "[regenotype]") {
    // Every read is at this site, of run 0, at one other site of run 0, and at one of run 1,
    // which says little. Leaving this site out, its strand is its run-0 part alone.
    RegenotypeParams params;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    const LambdaTable coherent = confident(ev);
    unordered_map<uint64_t, double> own;
    LambdaTable mixed;
    LambdaTable plain;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        const uint64_t key = ev.read_key[i];
        const double part = coherent.at(key).lambda;
        own[key] = 0.5 * part;
        ReadLambda rl;
        rl.phase_set = 1;
        rl.sites = 3;
        rl.runs[0] = 0;
        rl.run_sites[0] = 2;
        rl.run_lambda[0] = own[key] + part;
        rl.runs[1] = 1;
        rl.run_sites[1] = 1;
        rl.run_lambda[1] = -0.25 * part;
        rl.lambda = rl.run_lambda[0] + rl.run_lambda[1];
        mixed[key] = rl;
        ReadLambda one = rl;
        one.lambda = own[key] + part;
        one.runs[1] = NO_RUN;
        one.run_sites[1] = 0;
        one.run_lambda[1] = 0.0;
        one.sites = 2;
        plain[key] = one;
    }
    RegenotypeCounters counters;
    map<vector<int>, double> with_mixed = flat_gl();
    phase_aware_correction(ev, mixed, 1, 0, 0, own, StrandCalibration{1.0, 1.0}, params, with_mixed,
                           counters);
    REQUIRE(counters.sites_mixed_reads == 1);
    REQUIRE(counters.sites_split_runs == 0);
    map<vector<int>, double> without = flat_gl();
    phase_aware_correction(ev, plain, 1, 0, 0, own, StrandCalibration{1.0, 1.0}, params, without,
                           counters);
    REQUIRE(with_mixed.at({0, 1}) == Approx(without.at({0, 1})));
    REQUIRE(with_mixed.at({0, 1}) > flat_gl().at({0, 1}));
}

TEST_CASE("a lone read from another run adds nothing", "[regenotype]") {
    // Its order is averaged over, and with equal allele lengths the average of its two terms is
    // its unphased term, so the site scores as if the read had no strand.
    RegenotypeParams params;
    RegenotypeCounters counters;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable lone = in_run(confident(ev), 0);
    lone[ev.read_key[29]].runs[0] = 1;
    LambdaTable strandless = lone;
    strandless[ev.read_key[29]].lambda = 0.0;
    strandless[ev.read_key[29]].run_lambda[0] = 0.0;
    map<vector<int>, double> with_lone = flat_gl();
    phase_aware_correction(ev, lone, 1, 0, 0, {}, StrandCalibration{1.0, 1.0}, params, with_lone,
                           counters);
    REQUIRE(counters.sites_split_runs == 1);
    map<vector<int>, double> without = flat_gl();
    phase_aware_correction(ev, strandless, 1, 0, 0, {}, StrandCalibration{1.0, 1.0}, params,
                           without, counters);
    REQUIRE(with_lone.at({0, 1}) == Approx(without.at({0, 1})));
}

TEST_CASE("a read's run leaves out the site being corrected", "[regenotype]") {
    ReadLambda rl;
    rl.runs[0] = 2;
    rl.run_sites[0] = 3;
    rl.runs[1] = 5;
    rl.run_sites[1] = 1;
    REQUIRE(read_run(rl, NO_RUN) == MIXED_RUNS);
    REQUIRE(read_run(rl, 5) == 2);
    REQUIRE(read_run(rl, 2) == MIXED_RUNS);
    ReadLambda alone;
    alone.runs[0] = 2;
    alone.run_sites[0] = 1;
    REQUIRE(read_run(alone, NO_RUN) == 2);
    REQUIRE(read_run(alone, 2) == NO_RUN);
    REQUIRE(read_run(ReadLambda(), NO_RUN) == NO_RUN);
    rl.many_runs = true;
    REQUIRE(read_run(rl, 5) == MIXED_RUNS);
}

TEST_CASE("accumulate_lambda keeps each run's part of a read's log-odds", "[regenotype]") {
    // Read 7 is at three sites: two in run 0 and one in run 1.
    vector<PhaseSite> sites(3);
    for (size_t i = 0; i < sites.size(); ++i) {
        sites[i].record_key = i + 1;
        sites[i].phase_set = 1;
        sites[i].run = i < 2 ? 0 : 1;
        sites[i].read_key = {7};
        sites[i].q0 = {i == 2 ? 0.1f : 0.9f};
        sites[i].c = {0.9f};
    }
    LambdaTable lambda;
    RegenotypeCounters counters;
    accumulate_lambda(sites, {}, lambda, counters);
    const ReadLambda& rl = lambda.at(7);
    REQUIRE(rl.runs[0] == 0);
    REQUIRE(rl.run_sites[0] == 2);
    REQUIRE(rl.runs[1] == 1);
    REQUIRE(rl.run_sites[1] == 1);
    REQUIRE(rl.run_lambda[0] == Approx(2.0 * site_read_log_odds(0.9, 0.9)));
    REQUIRE(rl.run_lambda[1] == Approx(site_read_log_odds(0.1, 0.9)));
    REQUIRE(rl.lambda == Approx(rl.run_lambda[0] + rl.run_lambda[1]));
    REQUIRE(!rl.many_runs);
    // A fourth site in a third run is one run too many to keep apart.
    sites.push_back(sites[0]);
    sites.back().record_key = 4;
    sites.back().run = 2;
    LambdaTable more;
    accumulate_lambda(sites, {}, more, counters);
    REQUIRE(more.at(7).many_runs);
}

TEST_CASE("a nested chain counts the part of its own run from a read in two runs", "[regenotype]") {
    // Each read's run-0 part points it away from strand 0, where the chain is, and its larger
    // run-1 part points it toward strand 0. The chain's strand is in run 0's frame, so only the
    // run-0 part counts, and it makes the reads count less.
    RegenotypeParams params;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable lambda;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.phase_set = 1;
        rl.sites = 3;
        rl.runs[0] = 0;
        rl.run_sites[0] = 1;
        rl.run_lambda[0] = -4.0;
        rl.runs[1] = 1;
        rl.run_sites[1] = 2;
        rl.run_lambda[1] = 60.0;
        rl.lambda = 56.0;
        lambda[ev.read_key[i]] = rl;
    }
    const map<vector<int>, double> before = {{{0}, -20.0}, {{1}, -35.0}};
    RegenotypeCounters counters;
    map<vector<int>, double> own_run = before;
    haploid_inclusion_correction(ev, lambda, 1, 0, {}, StrandCalibration{1.0, 1.0}, +1, params,
                                 own_run, counters);
    REQUIRE(own_run != before);
    map<vector<int>, double> other_run = before;
    haploid_inclusion_correction(ev, lambda, 1, 1, {}, StrandCalibration{1.0, 1.0}, +1, params,
                                 other_run, counters);
    REQUIRE(other_run == before);
}

TEST_CASE("a nested chain leaves alone the reads whose strands come from another run",
          "[regenotype]") {
    // The chain's strand is named in the frame of run 0, its parent's. Reads whose strands come
    // from run 1 are related to that frame only by the panel.
    RegenotypeParams params;
    PhaseReadEvidence ev = evidence(30, 10, 10);
    LambdaTable lambda;
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        ReadLambda rl;
        rl.lambda = -60.0;
        rl.sites = 8;
        rl.phase_set = 1;
        lambda[ev.read_key[i]] = rl;
    }
    lambda = in_run(lambda, 1);
    const map<vector<int>, double> before = {{{0}, -20.0}, {{1}, -35.0}};
    RegenotypeCounters counters;
    map<vector<int>, double> other_run = before;
    haploid_inclusion_correction(ev, lambda, 1, 0, {}, StrandCalibration{1.0, 1.0}, +1, params,
                                 other_run, counters);
    REQUIRE(other_run == before);
    REQUIRE(counters.haploid_other_run == 1);
    map<vector<int>, double> same_run = before;
    haploid_inclusion_correction(ev, lambda, 1, 1, {}, StrandCalibration{1.0, 1.0}, +1, params,
                                 same_run, counters);
    map<vector<int>, double> unknown = before;
    haploid_inclusion_correction(ev, lambda, 1, NO_RUN, {}, StrandCalibration{1.0, 1.0}, +1,
                                 params, unknown, counters);
    REQUIRE(same_run != before);
    REQUIRE(unknown == same_run);
}

}
}
