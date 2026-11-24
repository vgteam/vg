#ifndef VG_READ_PHASING_HPP_INCLUDED
#define VG_READ_PHASING_HPP_INCLUDED

/** \file read_phasing.hpp
 *
 * Read phasing: re-decide the order of each heterozygous site's chosen allele pair from the reads
 * that span it and other heterozygous sites, within the phase sets the panel gave.
 *
 * Genotypes do not change. To flip a site is to swap its pair against the order the panel gave.
 * Each phase set is phased in stages:
 *
 *   1. Link each reliable site to the next, making the phase chain. Break it where a link is
 *      weak, and orient each piece from its first site.
 *   2. Decide each break from the links across it. Turn each run of pieces that reads join so
 *      that its first site the panel orders keeps the panel's order.
 *   3. Remove phase-chain sites whose reads disagree with the rest of the chain, and repeat
 *      stages 1 and 2.
 *   4. Hang every other site from its nearest phase-chain sites. A wrong decision there flips one
 *      site, where a wrong link in the chain would flip every site after it.
 *   5. Order by reads alone the sites the panel does not order that no read links to the chain.
 *
 * The method is described in artifacts/doc/read-likelihood-read-phasing.md, under "Read phasing".
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
 * A site's per-read evidence, kept from the direct pass until read phasing runs after the linkage
 * pass.
 *
 * It holds only what phasing needs, which is much less than AnchorSiteEvidence. A read is keyed
 * by a hash of its name, so paired mates, which share a name and lie on one haplotype, count as
 * one read. Two reads whose names collide are taken for one.
 */
struct PhaseReadEvidence {
    /// The key of each read.
    vector<uint64_t> read_key;
    /// The probability that each read is mismapped here.
    vector<float> mismap;
    /// rel(r, a), row major, reads x alleles, row-normalised into [0,1].
    vector<float> rel;
    /// The number of alleles, and of columns in `rel`.
    size_t n_alleles = 0;
    /// The full length of each allele, for the allele-length weights.
    vector<uint32_t> allele_length;
    /// The mean read length R the site's matrix used, for the allele-length weights.
    ///
    /// This is the rate window's mean; see AlleleReadLikelihoods::set_length_weights.
    float mean_read_length = 0.0f;
    /// Whether the allele-length weights are used, or flat weights instead.
    bool length_weighted = true;

    /// Get rel(read, allele).
    float rel_at(size_t read, size_t allele) const { return rel[read * n_alleles + allele]; }
    /// Get the number of reads.
    size_t num_reads() const { return read_key.size(); }
    /// Bytes held, for reporting.
    size_t bytes() const {
        return read_key.size() * sizeof(uint64_t) + mismap.size() * sizeof(float)
               + rel.size() * sizeof(float) + allele_length.size() * sizeof(uint32_t)
               + sizeof(PhaseReadEvidence);
    }
};

/**
 * Return the allele-length weights of a site's slots.
 *
 * Each weight is a slot's expected share of the site's reads, from the full length of the allele
 * it holds. An allele of length L yields a read of length R that overlaps the site from L + R - 1
 * start positions. `slot_allele` names the allele in each slot, so reordering the slots reorders
 * the weights. The weights are flat when `length_weighted` is false or the lengths are missing.
 */
vector<double> allele_length_weights(const vector<std::uint32_t>& allele_length, size_t n_alleles,
                                     float mean_read_length, bool length_weighted,
                                     const vector<int>& slot_allele);

/**
 * One heterozygous site's read evidence, reduced to its chosen pair.
 *
 * Each row is one read that fits at least one of the two alleles. The rows are sorted by read key,
 * one row per key (see `merge_mates`).
 */
struct PhaseSite {
    /// The key of the site's record, which names the site in the results.
    size_t record_key = 0;
    /// The site's phase set, as an id unique across contigs.
    size_t phase_set = 0;
    /// The site's position, which orders the sites of a phase set.
    size_t position = 0;
    /// The key of each read.
    vector<uint64_t> read_key;
    /// For each read, the probability that it carries the allele in slot 0 (strand 0), given that
    /// it came from one of the two chosen strands.
    vector<float> q0;
    /// For each read, the probability that it came from one of the two chosen strands, rather than
    /// being mismapped.
    vector<float> c;
    /// Mean of the rows' confidences (see `read_confidence`). Low where the reads cannot tell the
    /// two alleles apart. A site is reliable when this is at least `ReadPhasingParams::reliability`.
    double reliability = 0.0;
    /// Whether the panel orders the site's pair.
    ///
    /// It does not where no panel haplotype tells the strands apart (`PhaseCall::order_arbitrary`),
    /// as for a heterozygote whose ALT no haplotype carries. Such a site is ordered by its reads
    /// alone, and never lends the panel's order to another site.
    bool panel_ordered = true;
};

/// Return a read's confidence at a site.
///
/// This is the phred-scaled probability that its better allele is wrong,
/// `-10 log10(1 - max(q0, 1 - q0) * c)`, capped at 120.
double read_confidence(float q0, float c);

/// Sort a site's rows by read key and keep one row per key.
///
/// Paired mates share a key, and they come from one molecule and often read the same bases, so
/// their rows are not independent evidence. We keep the row with the higher confidence, then the
/// larger `q0`, then the larger `c`. The result does not depend on the order of the rows.
void merge_mates(PhaseSite& site);

/// Reduce a site's evidence to its chosen pair, with `allele0` in slot 0 and `allele1` in slot 1.
///
/// The rows are merged by `merge_mates`, and the reliability is computed over the merged rows. The
/// caller fills `record_key`, `phase_set` and `position`.
PhaseSite reduce_to_pair(const PhaseReadEvidence& evidence, size_t allele0, size_t allele1);

/// The settings of read phasing.
struct ReadPhasingParams {
    /// The reliability a site needs to be in the phase chain.
    ///
    /// Between two alleles of equal length, a read that fits one perfectly and the other not at
    /// all scores phred(e / (e + (1 - e) / 2)) for e the floor on the mismapping probability, the
    /// heterozygous score ceiling. Scores cluster just below it, so the threshold is sensitive to
    /// anything that moves them.
    double reliability = 9.5;
    /// The size of a link, in log10 units, below which the phase chain breaks.
    ///
    /// Every break is then decided by the relink, which weighs several links.
    double break_threshold = 20.0;
    /// The number of phase-chain sites on each side of a break whose links decide it.
    size_t relink = 10;
    /// The number of phase-chain sites a site off the chain is hung from.
    ///
    /// We use the nearest hang / 2 + 1 on each side, so up to 6 in all at the default of 4.
    size_t hang = 4;
    /// The weight, in log10 units, of a vote for the panel's order when hanging a site.
    ///
    /// At 0, a site that no read links to the chain keeps the panel's order.
    double panel_weight = 3.0;
    /// The limit on the size of one link, or 0 for none.
    ///
    /// Reads are treated as independent, so a link from many reads can be far larger than its real
    /// reliability.
    double cap = 0.0;

    /// The coherence a site needs to stay in the phase chain, in [0,1], or 0 for no coherence stage.
    ///
    /// A site's coherence is the fraction of its reads whose allele agrees with the strand that
    /// the read's other phase-chain sites put it on. A read agrees when its q0 is on the same side
    /// of 1/2 as that strand. Reliability says whether a site's reads can tell its alleles apart,
    /// and coherence whether they agree with the rest of the chain.
    double coherence_min = 0.70;
    /// The most rounds of the coherence stage that may remove sites.
    ///
    /// Each round measures coherence on the phase chain the previous round produced. More rounds
    /// remove more sites, so the remaining links span further, are weaker, and break more often.
    size_t coherence_rounds = 2;
};

/// A report of what read phasing did.
struct ReadPhasingCounters {
    /// Sites phased.
    size_t sites = 0;
    /// Reliable sites.
    size_t reliable = 0;
    /// Phase sets, each with one phase chain.
    size_t chains = 0;
    /// Breaks in the phase chains.
    size_t breaks = 0;
    /// Breaks that no read links across.
    size_t breaks_no_reads = 0;
    /// Sites hung from a phase chain.
    size_t hung = 0;
    /// Hung sites that no read ordered.
    size_t hung_no_reads = 0;
    /// Sites flipped.
    size_t flipped = 0;
    /// Nested sites whose `nested_strand` was inverted because their parent's pair was swapped.
    size_t strands_rederived = 0;
    /// Sites removed from the phase chain for low coherence.
    size_t demoted_incoherent = 0;
    /// The most rounds of the coherence stage that removed sites from one phase chain.
    size_t coherence_rounds_run = 0;
    /// Phase chains still losing sites in the last round allowed.
    size_t coherence_unconverged = 0;
};

/// Return the size of the link between two sites.
///
/// This is the log10 odds, cis against trans, over the reads the sites share. Positive means the
/// reads agree with the sites' current slot orders.
double phase_link(const PhaseSite& a, const PhaseSite& b, double cap);

/**
 * Decide which sites to flip.
 *
 * Phase is comparable only inside a phase set, so we sort `sites` by `phase_set`, then
 * `position`, then `record_key`, in place. Each site's rows are merged by `merge_mates`, in place,
 * so a site built by hand is treated as `reduce_to_pair` would build it.
 *
 * With no read evidence, no site is flipped.
 *
 * @returns The record keys of the sites to flip.
 *
 * @param read_ordered If not null, receives the record keys of the sites the panel does not order
 * whose reads ordered them.
 */
unordered_set<size_t> read_phase_flips(vector<PhaseSite>& sites, const ReadPhasingParams& params,
                                       ReadPhasingCounters& counters,
                                       unordered_set<size_t>* read_ordered = nullptr);

}

#endif
