/// \file unittest/read_phasing.cpp
///
/// Unit tests for the three-stage read-backed phase decision.
///
/// Everything here is parity arithmetic over a chain, which is the class of error that produces a
/// plausible-looking phase and a wrong haplotype: a sign slip downstream of one link flips every
/// site after it and nothing in the VCF looks unusual. The chr20 switch error catches it on real
/// data; these pin the convention so that when that number moves it means the data surprised us
/// rather than that the algebra was never right.
///

#include <vector>

#include "catch.hpp"
#include "read_phasing.hpp"

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
        s.p.push_back(0.95f);
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
