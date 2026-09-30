/// \file unittest/regenotype.cpp
///
/// Unit tests for the phase-aware genotype correction.
///
/// Neither `--regeno-temper 0` nor `--regeno-passes 1` exercises `Lambda`: at tau = 0 the
/// weighting is the identity whatever `Lambda` holds, and with one pass nothing is recorded again.
/// Taking the larger correction over the two orders also hides a `Lambda` whose sign is wrong at
/// every site. So `Lambda` wired backwards at some sites would look like the shuffled control, and
/// these tests are what check its sign.

#include <cmath>
#include <map>
#include <vector>

#include "catch.hpp"
#include "anchor.hpp"
#include "regenotype.hpp"

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
    // The gate the whole arm rests on, and it holds by construction rather than by tolerance:
    // at tau = 0 the tilted weights ARE the slot weights, so the corrected and uncorrected
    // mixtures are the same arithmetic on the same values.
    RegenotypeParams params;
    RegenotypeCounters counters;
    for (uint32_t len1 : {10u, 400u, 12000u}) {
        PhaseReadEvidence ev = evidence(30, 10, len1);
        LambdaTable lambda = confident(ev);
        map<vector<int>, double> gl = flat_gl();
        const map<vector<int>, double> before = gl;
        phase_aware_correction(ev, lambda, {}, 0.0, 1.0, params, gl, counters);
        for (const auto& kv : before) {
            REQUIRE(gl.at(kv.first) == kv.second);
        }
    }
    // And the two orders scored EQUAL, not merely close. If they can differ by a last bit at
    // temper 0 then the correction is not exactly zero, it is zero-ish -- which is a different
    // and much weaker claim than the one the byte-identity gate is asked to carry.
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
        phase_aware_correction(ev, lambda, {}, tau, 1.0, params, gl, counters);
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
    phase_aware_correction(ev, lambda, {}, 2.0, 1.0, params, gl, counters);
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
    // flip set describes the panel's order, not the settled one, and its contribution must enter
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
    // rows under one key. A fragment lies on one strand, so it counts once. `phase_link`'s sorted
    // merge pairs them one to one; a running sum would count them twice.
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
    phase_aware_correction(ev, coherent, {}, 1.0, 1.0, params, gl_coherent, counters);

    // Same magnitudes, strand assignment uncorrelated with the allele carried.
    LambdaTable incoherent = confident(ev);
    for (size_t i = 0; i < ev.num_reads(); ++i) {
        incoherent[ev.read_key[i]].lambda = ((i / 2) % 2 == 0) ? 4.0 : -4.0;
    }
    map<vector<int>, double> gl_incoherent = flat_gl();
    phase_aware_correction(ev, incoherent, {}, 1.0, 1.0, params, gl_incoherent, counters);

    const double het_gap_coherent = gl_coherent.at({0, 1}) - gl_coherent.at({0, 0});
    const double het_gap_incoherent = gl_incoherent.at({0, 1}) - gl_incoherent.at({0, 0});
    const double het_gap_before = -20.0 - -50.0;
    REQUIRE(het_gap_coherent > het_gap_before);
    REQUIRE(het_gap_incoherent < het_gap_coherent);
}

TEST_CASE("the per-read escape is exactly nothing at temper 0, whatever the ceiling",
          "[regenotype]") {
    // The escape has to leave the byte-identity gate alone, and it does so for a reason that does
    // not depend on the fitted value: at tau = 0 the sigmoid is 1/2 for any ceiling, so the
    // probability is 1/2 and its logit is 0.
    for (double c : {1.0, 0.99, 0.95, 0.5, 0.01}) {
        for (double lam : {0.0, 17.3, -1215.0, 1e9, -1e-9}) {
            REQUIRE(calibrated_log_odds(lam, 0.0, c) == 0.0);
        }
    }
}

TEST_CASE("the per-read escape caps how confident a read may be", "[regenotype]") {
    // The point of the ceiling. At c = 1 a read at |Lambda| = 1215 and the fitted temper asserts
    // a strand at odds of e^60, which the calibration says is wrong -- observed agreement there
    // is 0.973, not 1. The escape holds it to something the data supports.
    const double tau = 0.0498438;
    const double wild = calibrated_log_odds(1215.0, tau, 1.0);
    const double held = calibrated_log_odds(1215.0, tau, 0.95);
    REQUIRE(wild > 20.0);
    REQUIRE(held < 5.0);
    REQUIRE(std::isfinite(wild));   // clamped off 1, or the logit is infinite
    // Monotone in the evidence, and still ordered the same way.
    REQUIRE(calibrated_log_odds(102.0, tau, 0.95) < held);
    REQUIRE(calibrated_log_odds(-1215.0, tau, 0.95) == Approx(-held).epsilon(1e-12));
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
        haploid_inclusion_correction(ev, lambda, {}, 0.0, 0.95, sign, params, gl, counters);
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
    // The likelihoods must be the read term, as the sweep's are: the correction subtracts the read
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
    haploid_inclusion_correction(ev, lambda, {}, 1.0, 1.0, +1, params, gl, counters);
    // Both alleles rise, since reads that do not belong here no longer count against them. The
    // correction is added to the sweep's likelihood, so the result is not bounded above by zero.
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
    haploid_inclusion_correction(ev, here, {}, 1.0, 1.0, +1, params, kept, counters);
    REQUIRE(kept.at({0}) == kept_before.at({0}));
    REQUIRE(kept.at({1}) == kept_before.at({1}));
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
    phase_aware_correction(ev, normal, {}, 1.0, 1.0, params, a, counters);
    phase_aware_correction(ev, flipped, {}, 1.0, 1.0, params, b, counters);
    REQUIRE(a.at({0, 1}) == Approx(b.at({0, 1})).epsilon(1e-12));
}

}
}
