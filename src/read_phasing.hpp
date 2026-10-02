#ifndef VG_READ_PHASING_HPP_INCLUDED
#define VG_READ_PHASING_HPP_INCLUDED

/** \file read_phasing.hpp
 *
 * Read-backed phasing (--read-phasing): re-decide the order of each heterozygous site's settled
 * allele pair from the reads that span it and other heterozygous sites, within the phase sets the
 * panel gave. No genotype changes.
 *
 * The orders are decided in stages:
 *
 *   1. Chain the reliable sites of a phase set, each linked to the next reliable site, so that
 *      unreliable sites between them are stepped over. Break the chain where a link is weak.
 *   2. Relink each break from the reliable sites on either side of it.
 *   3. Remove chain sites whose reads disagree with the strand the chain gives them (low
 *      coherence), and repeat stages 1 and 2.
 *   4. Hang every other site from its nearest chain sites. A wrong decision there flips one site,
 *      where a wrong link in the chain would flip every site after it.
 *
 * The method is described in doc/read-likelihood-genotyping.md, under "From the reads".
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
 * A site's per-read evidence, kept from the sweep until read phasing runs after the barrier.
 *
 * It holds only what phasing needs, which is much less than AnchorSiteEvidence. A read is keyed
 * by a hash of its name, since paired mates share a name and lie on one haplotype, so they count
 * as one read. Two reads whose names collide are taken for one.
 */
struct PhaseReadEvidence {
    vector<uint64_t> read_key;
    vector<float> mismap;
    /// rel(r, a), row major, reads x alleles, row-normalised into [0,1].
    vector<float> rel;
    size_t n_alleles = 0;
    /// The alleles' full lengths, and the mean read length R the site's matrix used (the rate
    /// window's mean; see AlleleReadLikelihoods::set_length_weights), for the allele-length
    /// weights.
    vector<uint32_t> allele_length;
    float mean_read_length = 0.0f;
    bool length_weighted = true;

    float rel_at(size_t read, size_t allele) const { return rel[read * n_alleles + allele]; }
    size_t num_reads() const { return read_key.size(); }
    /// Bytes held, for reporting.
    size_t bytes() const {
        return read_key.size() * sizeof(uint64_t) + mismap.size() * sizeof(float)
               + rel.size() * sizeof(float) + allele_length.size() * sizeof(uint32_t)
               + sizeof(PhaseReadEvidence);
    }
};

/// The allele-length weights of a site's slots: each slot's expected share of the site's reads,
/// from the full length of the allele it holds. An allele of length L yields a read of length R
/// that overlaps the site from L + R - 1 start positions. `slot_allele` names the allele in each
/// slot, so reordering the slots reorders the weights. The weights are flat when
/// `length_weighted` is false or the lengths are missing.
vector<double> allele_length_weights(const vector<std::uint32_t>& allele_length, size_t n_alleles,
                                     float mean_read_length, bool length_weighted,
                                     const vector<int>& slot_allele);

/**
 * One heterozygous site's read evidence, reduced to its settled pair.
 *
 * Each row is one read. `q0` is the probability that the read carries the allele in slot 0
 * (strand 0), given that it came from one of the two settled strands, and `c` is the probability
 * that it came from one of them, rather than being mismapped. Only reads that fit one of the two
 * alleles at all are kept. The rows are sorted by read key, one row per key (see `merge_mates`).
 */
struct PhaseSite {
    size_t record_key = 0;
    /// The site's phase set, as an id unique across contigs.
    size_t phase_set = 0;
    size_t position = 0;
    vector<uint64_t> read_key;
    vector<float> q0;
    vector<float> c;
    /// Mean of the rows' confidences (see `read_confidence`). Low where the reads cannot tell the
    /// two alleles apart. A site is reliable when this is at least `ReadPhasingParams::reliability`.
    double reliability = 0.0;
};

/// A read's confidence at a site: the phred-scaled probability that its better allele is wrong,
/// `-10 log10(1 - max(q0, 1 - q0) * c)`, capped at 120.
double read_confidence(float q0, float c);

/// Sort a site's rows by read key and keep one row per key.
///
/// Paired mates share a key, so where both reach the site they are one read. They come from one
/// molecule and often read the same bases, so their rows are not independent evidence, and only
/// the row with the higher confidence is kept; between equal confidences, the larger `q0`, then
/// the larger `c`. The result does not depend on the order of the rows.
void merge_mates(PhaseSite& site);

/// Reduce a site's evidence to its settled pair, `allele0` in slot 0 and `allele1` in slot 1: the
/// rows, merged by `merge_mates`, and the site's reliability over the merged rows. The caller fills
/// `record_key`, `phase_set` and `position`.
PhaseSite reduce_to_pair(const PhaseReadEvidence& evidence, size_t allele0, size_t allele1);

struct ReadPhasingParams {
    /// A site below this reliability cannot be in the chain (--phase-min-q).
    ///
    /// Between two alleles of equal length, a read that fits one perfectly and the other not at
    /// all scores phred(e / (e + (1 - e) / 2)) for e = --mismap-min, the heterozygous score
    /// ceiling. Scores cluster just below it, so the threshold is sensitive to anything that
    /// moves them.
    double reliability = 9.5;
    /// Break the chain where the size of a link is below this, in log10 units (--phase-break).
    /// Every break is then decided by the relink, which weighs several links.
    double break_threshold = 20.0;
    /// Chain sites on each side of a break whose links decide it (--phase-relink).
    size_t relink = 10;
    /// Which chain sites an unreliable site is hung from (--phase-hang): the nearest
    /// hang / 2 + 1 on each side, so up to 6 in all at the default of 4.
    size_t hang = 4;
    /// Weight, in log10 units, of a vote for the panel's order when hanging a site
    /// (--phase-prior). At 0, a site that no read links to the chain keeps the panel's order.
    double panel_weight = 3.0;
    /// Limit on the size of one link, 0 for none (--phase-cap). Reads are treated as independent,
    /// so a link from many reads can be far larger than its real reliability.
    double cap = 0.0;

    /// Minimum coherence for a site to stay in the chain, in [0,1]; 0 turns the coherence
    /// step off (--phase-coherence).
    ///
    /// A site's coherence is the fraction of its reads whose allele agrees with the strand that
    /// the read's other chain sites put it on, leaving the site itself out; a read agrees when
    /// its q0 is on the same side of 1/2 as that strand. Reliability asks whether a site's reads
    /// can tell its alleles apart; coherence asks whether they agree with the rest of the chain.
    /// Sites below the minimum, among those with enough such reads, are marked unreliable, and
    /// the chain and relink steps run again.
    double coherence_min = 0.70;
    /// The most rounds of the coherence step (--phase-coh-rounds). Each round measures coherence on the
    /// chain the previous round produced. More rounds remove more sites, so the remaining links
    /// span further, are weaker, and break more often.
    size_t coherence_rounds = 2;

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
    size_t strands_rederived = 0;
    /// Sites removed from the chain for low coherence.
    size_t demoted_incoherent = 0;
    /// The most rounds of the coherence step that removed sites from one chain, and the chains
    /// still losing sites when the round limit was reached.
    size_t coherence_rounds_run = 0;
    size_t coherence_unconverged = 0;
};

/// log10 odds, cis against trans, over the reads two sites share. Positive means the reads agree
/// with the sites' current slot order.
double phase_link(const PhaseSite& a, const PhaseSite& b, double cap);

/// Decide every site's order. `sites` may arrive in any order; they are grouped by `phase_set`
/// and sorted by `position` internally, since phase is comparable only inside a phase set. Each
/// site's rows are merged by `merge_mates`, in place, so a site built by hand is treated as
/// `reduce_to_pair` would build it.
///
/// Returns the record keys whose settled pair should be swapped. A chain's first site keeps the
/// panel's order, so with no read evidence nothing is swapped. Breaks are relinked from left to
/// right, each piece taking its orientation from the piece before it as that piece now stands. A
/// hung site whose votes sum to exactly zero keeps the panel's order.
unordered_set<size_t> read_phase_flips(vector<PhaseSite>& sites, const ReadPhasingParams& params,
                                       ReadPhasingCounters& counters);

}

#endif
