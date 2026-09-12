#ifndef VG_READ_PHASING_HPP_INCLUDED
#define VG_READ_PHASING_HPP_INCLUDED

/** \file read_phasing.hpp
 *
 * Decide each heterozygous site's phase from the reads that cross it.
 *
 * The phase `vg call --phased` ships comes from the haplotype panel: the linkage layer's Viterbi
 * path names a pair of panel haplotypes and the order of that pair is the phase. The reads are never
 * asked. On chr20 ONT 44x that costs a 3.79% switch error while the reads in the same run answer at
 * 0.09% across the sites they can speak for.
 *
 * Three stages, and adjacency is deliberately NOT the axis:
 *
 *   1. link *consecutive reliable* sites, stepping over unreliable ones so a homopolymer never sits
 *      inside the chain it would break. Break the chain where the link is weak.
 *   2. relink each break over a few reliable sites either side.
 *   3. hang the unreliable sites off the settled chain, where a misplacement is a FLIP -- two
 *      junction errors that do not propagate -- rather than a switch that propagates forever.
 *
 * MEASURED, and it is the reason for the shape. At the junctions stage 1 breaks, the adjacent link
 * is 15.06% wrong; the same junction decided by the nearest *reliable* site either side is 0.97%
 * wrong; but the nearest site regardless of quality is 14.93% wrong -- indistinguishable from
 * adjacent. Reaching further buys nothing. Stepping over bad sites is the entire effect.
 *
 * The evidence is not circular with the panel: a read's allele responsibility depends on the settled
 * *genotype*, which is phase-independent. The phase only relabels the slots.
 */

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

namespace vg {

using std::size_t;
using std::uint64_t;
using std::unordered_set;
using std::vector;

/**
 * A site's per-read evidence, retained from the read sweep to the barrier.
 *
 * Deliberately NOT `AnchorSiteEvidence`, which the first implementation reused. That struct carries
 * a read NAME -- a 36-character ONT UUID, so a heap allocation per read per site -- plus two
 * resolved pins, together about 120 bytes a read against the 24 this needs. Retaining it for
 * phasing cost **+2.8 GB of peak RSS on chr6** (7.1 -> 9.9 GB).
 *
 * The key is a hash of the read name, which is the right identity rather than a compromise: paired
 * mates share a name and a fragment lies on one haplotype, so one name is one observation. A
 * collision mislabels one read at one pair and degrades gracefully, which a lock on a global intern
 * table in the sweep's hot path would not.
 */
struct PhaseReadEvidence {
    vector<uint64_t> read_key;
    vector<float> mismap;
    /// rel(r, a), row major, reads x alleles, row-normalised into [0,1].
    vector<float> rel;
    size_t n_alleles = 0;
    /// The alleles' spelled lengths and the site's mean read length, for the mixture weights.
    vector<uint32_t> allele_length;
    float mean_read_length = 0.0f;
    bool length_weighted = true;

    float rel_at(size_t read, size_t allele) const { return rel[read * n_alleles + allele]; }
    size_t num_reads() const { return read_key.size(); }
    /// Retained bytes, reported rather than estimated.
    size_t bytes() const {
        return read_key.size() * sizeof(uint64_t) + mismap.size() * sizeof(float)
               + rel.size() * sizeof(float) + allele_length.size() * sizeof(uint32_t)
               + sizeof(PhaseReadEvidence);
    }
};

/**
 * One heterozygous site's read evidence, already reduced to the settled pair.
 *
 * Per read, two numbers are enough. `q0` is the read's probability of carrying slot 0's allele
 * *given* that it came from one of the two settled haplotypes, and `p` is the probability that it
 * came from one of them at all -- one minus the mismapping escape the genotype model already
 * clamps. Nothing else about the read survives, and nothing else is needed.
 */
struct PhaseSite {
    size_t record_key = 0;
    size_t phase_set = 0;
    size_t position = 0;
    vector<uint64_t> read_key;
    vector<float> q0;
    vector<float> p;
    /// Mean over reads of the phred complement of the winner's share -- the same quantity the anchor
    /// file writes per placement. Low exactly where the reads cannot tell which allele they carry,
    /// which on ONT means a 1 bp indel: 99.67% of those sites fall below the default threshold
    /// against 2.14% of SNVs.
    double reliability = 0.0;
};

struct ReadPhasingParams {
    /// A site below this is not allowed to carry a link. Fitted on chr20; the distribution is tight
    /// (median 10.09, 25th percentile 9.95) so this is sensitive and wants re-fitting whenever the
    /// per-read scores move. `--mismap-min` moves them directly: the score's ceiling is phred of
    /// it, 13.01 at the preset's 0.05 against 16.99 at the 0.02 default.
    double reliability = 9.5;
    /// Break the chain below this many log10 units of evidence.
    double break_threshold = 10.0;
    /// Reliable sites either side of a break to relink over. 3 is enough; 8 is a wash and 15 hurt.
    size_t relink = 3;
    /// Decided neighbours to hang an unreliable site from.
    size_t hang = 4;
    /// Weight of the panel's own answer when hanging a site. The panel is a genuinely long-range
    /// prior and this is the one place it still earns its keep: it moved 0.501% to 0.477%.
    double panel_weight = 3.0;
    /// Clamp on one pair's contribution, 0 to disable. 40 independent reads give 10^40 while the
    /// measured read-only error is 10^-2.4, so the independence assumption is wrong by orders of
    /// magnitude and one chimeric link in a segmental duplication would otherwise be unoverridable.
    double cap = 0.0;
};

struct ReadPhasingCounters {
    size_t sites = 0;
    size_t reliable = 0;
    size_t chains = 0;
    size_t breaks = 0;
    size_t breaks_no_reads = 0;
    size_t hung = 0;
    size_t hung_no_reads = 0;
    size_t flipped = 0;
    /// Nested sites whose `nested_strand` was inverted because their parent's pair was swapped.
    /// Zero when nothing was re-phased; a re-phase that moves no strand at all under `-A` means the
    /// cascade is not reaching the tree.
    size_t strands_rederived = 0;
};

/// log10 odds, cis against trans, over the reads two sites share. Positive means the reads agree
/// with the sites' current slot order.
double phase_link(const PhaseSite& a, const PhaseSite& b, double cap);

/// Decide every site's orientation. `sites` may arrive in any order; it is grouped by `phase_set`
/// and sorted by `position` internally, because phase is only comparable inside a block.
///
/// Returns the record keys whose settled pair should be swapped. A chain's first site is never
/// swapped, which keeps the panel's global frame -- arbitrary either way -- and makes "no read
/// evidence anywhere" return the empty set rather than an inverted genome.
unordered_set<size_t> read_phase_flips(vector<PhaseSite>& sites, const ReadPhasingParams& params,
                                       ReadPhasingCounters& counters);

}

#endif
