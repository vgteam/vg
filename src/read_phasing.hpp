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
    /// per-read scores move. `--mismap-min` moves them directly.
    ///
    /// **The ceiling that applies here is the HETEROZYGOUS one, and it is not phred(--mismap-min).**
    /// A PhaseSite is only ever built for a diploid heterozygote (`FlowCaller::apply_read_phasing`
    /// skips anything with `trav_first == trav_second`), and a balanced het splits the slot weight
    /// in half, so a perfectly discriminating read reaches
    ///
    ///     phred( e / (e + (1 - e) / 2) )  =  10.21 at e = 0.05,  14.07 at e = 0.02
    ///
    /// against phred(e) = 13.01 and 16.99, which is the HOMOZYGOUS ceiling -- there the half weight
    /// is deliberately not applied, so it is the figure `--anchors-min-q` quotes and it does not
    /// apply to this gate. Measured on chr20 ONT: het sites run a median 8.98, p99 10.21 and a max
    /// of 10.36, while single-slot sites sit flat on 13.01.
    ///
    /// Two consequences. The operating range is about [2.8, 10.4], so this default sits under two
    /// phred from a hard cap on a nearly saturated statistic -- the reliable/unreliable split is a
    /// knife edge, not a comfortable classification. And any value above the het ceiling leaves
    /// `rel` empty in every block, which skips all three stages silently: `main_call` refuses such
    /// a value rather than letting a run report "0 reliable" and no phasing.
    ///
    /// That warning came true. `--realign` moved the distribution's median to 8.98 and left this
    /// 9.5 above almost all of it: 5.9% of heterozygous sites stayed eligible against 77.4%, and
    /// chr20 switch error went 0.3545% -> 0.5794%. `main_call` therefore overrides this to 8.5
    /// when the exact walk is in use; see the comment there and docs/phase-min-q-refit.md. This
    /// default is the GREEDY walk's value and is still right for it.
    double reliability = 9.5;
    /// Break the chain below this many log10 units of evidence.
    ///
    /// Raised from 10 after two things were established. First, a chain break does NOT fragment the
    /// output: 707 breaks and 7,565 breaks both give ONE `PS` block for a whole contig, because
    /// stage 2 relinks every break. So a break is not a cost, it is a REROUTE -- the junction stops
    /// being decided by stage 1's single-link sign cascade and starts being decided by stage 2's
    /// nine-pair magnitude-weighted relink. Second, the relink is the better rule, which is the
    /// same lesson `--phase-coherence` teaches by moving sites to stage 3's hanging rule.
    ///
    /// So breaking MORE often is breaking more often into a better decision procedure. Measured:
    /// chr20 33 -> 31 switches and chr6 31 -> 30 for this parameter alone, and 33 -> 27 / 31 -> 28
    /// together with the two below. Going the other way is worse -- `--phase-break 5` gives chr6 33
    /// against 31, which is the mechanism confirming itself.
    double break_threshold = 20.0;
    /// Reliable sites either side of a break to relink over.
    ///
    /// The old note here read "3 is enough; 8 is a wash and 15 hurt". That was measured BEFORE
    /// `--phase-coherence`, which removes from the backbone exactly the sites a wider window used
    /// to reach into, and before `--phase-break` was raised, which is what makes the relink carry
    /// most of the chain instead of a twentieth of it.
    ///
    /// Held to a lower standard of evidence than the other two, and that should be recorded: on
    /// SWITCHES it has no consistent direction (chr20 33 -> 29 alone but no change in combination;
    /// chr6 nothing alone but 30 -> 28 in combination). What argues for it is F1, on a denominator
    /// of 269,660 true positives rather than thirty switch events -- it is the only arm that moves
    /// F1 at all, and it moves it UP on chr6 by 0.00031, with TP +80 and FP -97, both directions
    /// right. One observation on one contig; a chr6 repeat would settle whether it is real or the
    /// genotype perturbation any phasing change causes through re-genotyping.
    size_t relink = 10;
    /// Decided neighbours to hang an unreliable site from.
    size_t hang = 4;
    /// Weight of the panel's own answer when hanging a site. The panel is a genuinely long-range
    /// prior and this is the one place it still earns its keep: it moved 0.501% to 0.477%.
    double panel_weight = 3.0;
    /// Clamp on one pair's contribution, 0 to disable. 40 independent reads give 10^40 while the
    /// measured read-only error is 10^-2.4, so the independence assumption is wrong by orders of
    /// magnitude and one chimeric link in a segmental duplication would otherwise be unoverridable.
    double cap = 0.0;

    /// Minimum aggregate log10 gain before a junction is treated as a switch and everything
    /// downstream of it flipped. 0 disables the pass. PROTOTYPE, off by default.
    ///
    /// Stage 1 decides a junction from the reads the two ADJACENT sites share, and uses only the
    /// sign. This pass instead asks, for every junction, what flipping all downstream orientations
    /// would do to the total read likelihood, using each read's ENTIRE span rather than one pair:
    ///
    ///     S(j) = sum over reads spanning j of sum over that read's sites s >= j of (b_s - a_s)
    ///
    /// where `a` is log10 of the read's per-site term under the haplotype it is assigned and `b`
    /// the term under the other. S(j) > 0 means the flip raises sum_r max_H L(H, r), so greedily
    /// flipping at the argmax junction is hill-climbing on that objective and terminates. This is
    /// the transitive constraint `phase_link` discards: a read spanning three sites constrains the
    /// 1->3 relation, and nothing in stages 1-3 ever reads it.
    double changepoint_min = 0.0;
    /// Cap on greedy flips per chain, so a pathological locus cannot spin.
    size_t changepoint_rounds = 200;

    /// Minimum PHASE COHERENCE for a site to keep carrying a link, in [0,1]. 0 disables.
    ///
    /// GATED ON BOTH CONTIGS AND ON BOTH METRICS, which no previous switch intervention managed:
    ///
    ///                     chr20 switches   chr6 switches   chr20 ALL F1   chr6 ALL F1
    ///     off                   51              62           0.95848        0.96487
    ///     0.70                  33              31           0.95900        0.96517
    ///
    /// chr6 is the hold-out and was not used to choose 0.70; chr20's own all-switch rate picks it
    /// (0.3142% against 0.3261% at 0.60 and 0.3312% at 0.80) and chr6 agrees independently. The
    /// assessed denominator RISES on both -- 58,799 -> 58,880 and 160,054 -> 160,137 -- so this is
    /// not whatshap's shrinking-denominator artefact, and F1 moves with TP up AND FP down on both,
    /// so it is not a precision/recall trade either.
    ///
    /// Coherence asks a different question from `reliability`. Reliability asks whether a site's
    /// reads can tell its two alleles apart; coherence asks whether those reads agree with the
    /// haplotype their OTHER sites imply. A site can be perfectly discriminable and completely
    /// phase-incoherent -- which is exactly the site a link must not be built from, and exactly the
    /// site reliability cannot see. Measured on chr20 ONT against switch positions, ranking sites
    /// worst-first: coherence enriches 9.3x in its worst 0.1% and 5.1x in its worst 1%, while
    /// reliability manages 2.2x and 1.1x and is BELOW the base rate at 5%. Pearson r between them
    /// is 0.437.
    ///
    /// Held out, not circular: a read's haplotype is recomputed for each site with that site's own
    /// term removed, so a site never votes on itself. Sites below the bar are demoted to unreliable
    /// and the cascade re-run, which is one extra pass.
    double coherence_min = 0.70;
    /// Demotion rounds before giving up on a coherent chain. 1 is the one-shot heuristic: demote
    /// once against the FIRST phasing and re-derive. Higher iterates -- re-measure coherence on the
    /// chain that re-derivation produced, demote again, repeat -- so the fixed point is a reliable
    /// chain every one of whose sites is coherent WITH THAT CHAIN, which the one-shot version does
    /// not guarantee. The risk is fragmentation: each round removes sites, the surviving links span
    /// further, `phase_link` falls off with distance, and more of them drop under --phase-break.
    ///
    /// 2 rather than 1: chr20 33 -> 31 switches and chr6 31 -> 30. Not the fixed point -- that is
    /// reached in three rounds and scores what one round scores, for 41% more chain breaks -- so
    /// this is a measured waypoint and not a principled limit, which is why it is capped here
    /// rather than run to convergence.
    size_t coherence_rounds = 2;
    /// Reads a site needs before low coherence may demote it.
    ///
    /// The instinct behind a high value is that coherence on few reads is noisy. The counter-case
    /// is stronger: a site with three reads that disagree contributes almost nothing to the
    /// backbone AND is unreliable, so protecting it is backwards -- it is a prime candidate, not a
    /// borderline one. Swept, because the argument does not settle it.
    size_t coherence_min_reads = 10;

    /// Backbone neighbours per side for the local-search refinement, 0 to disable.
    ///
    /// Stage 1 chains ADJACENT reliable sites and reads only the SIGN of each link, so one wrong
    /// sign inverts every site to the end of the segment -- mean 55.5 of them. Its magnitude is
    /// used once, for the break test, and a lone link has no competing evidence to weight against,
    /// which is why --phase-break was measured at 20/40/80 and changes nothing.
    ///
    /// This joins each backbone site to its K nearest backbone neighbours on EACH side and chooses
    /// orientations by local search on
    ///
    ///     sum over edges of |d_ij| * (+1 if the pair agrees with sign d_ij, else -1)
    ///
    /// flipping any site whose incident disagreeing weight exceeds its agreeing weight. Every flip
    /// strictly increases a bounded objective, so it terminates. It is stage 3's rule -- several
    /// neighbours, weighted by evidence -- applied to the backbone itself rather than only to the
    /// sites that were kept out of it.
    ///
    /// Single-site moves cannot cross a switch, which needs a whole suffix inverted; --phase-cp
    /// supplies that move, and its firing on 3 junctions in all of chr20 says suffix moves are the
    /// safety net and these are where the gain should be.
    size_t backbone = 0;

    /// Previous decided sites to consult when adding the NEXT site to the chain, 0 for the plain
    /// adjacent-pair cascade.
    ///
    /// Stage 1 decides each site from the single link to its predecessor, and from that link's SIGN
    /// alone, so one bad link inverts every site to the end of the segment -- mean 55.5 of them.
    /// With lookback the new site is decided by a weighted vote over the previous K sites already
    /// in the chain: sum over j of |d(m+1-j, m+1)| * (+1 if that site's settled orientation implies
    /// o = 0 for the new one, else -1). A single bad link is then outvoted by its neighbours instead
    /// of being obeyed.
    ///
    /// This is stage 3's rule used PROSPECTIVELY, during construction, rather than retrospectively
    /// on the sites the chain would not take. That distinction is why it is worth trying even
    /// though --phase-backbone was negative: local search refines the cascade's own answer and so
    /// settles into the nearest optimum, which is exactly what it found. A different construction
    /// lands in a different basin.
    ///
    /// Lookback never crosses a break: a new segment restarts at o = 0 with nothing behind it.
    size_t lookback = 0;

    /// Minimum weighted TRIANGLE CONSISTENCY for a site to enter the backbone at all, 0 to disable.
    ///
    /// Every other gate here is post-hoc: build the chain over all reliable sites, then judge each
    /// site against the chain it is already part of. That is circular, and it is why iterating the
    /// coherence demotion converges to no improvement -- the contaminated backbone is the yardstick.
    ///
    /// A triangle is frame-free. For sites i, j, k the reads imply sign(d_ij)*sign(d_jk)*sign(d_ik)
    /// > 0: going round a loop must flip an even number of times, whatever orientation anything is
    /// eventually given. So a site's consistency with its neighbourhood is measurable BEFORE any
    /// phase exists. Score each site by the share of its triangles that close, weighted by the
    /// weakest link in each so that a triangle resting on a marginal link cannot dominate, and admit
    /// only sites above the bar. The rest are hung by stage 3, exactly as a demoted site is.
    ///
    /// The object being sought is the largest set of sites that are strongly phase-consistent with
    /// ONE ANOTHER, rather than the largest set that individually separate their own alleles, which
    /// is all --phase-min-q can ask.
    double triangle_min = 0.0;
    /// Neighbours each side to draw triangles from.
    size_t triangle_k = 4;
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
    /// Junctions the changepoint pass flipped, and the aggregate log10 gain it claimed for them.
    size_t changepoints = 0;
    double changepoint_gain = 0.0;
    /// Junctions that cleared the threshold but were refused because the pass hit its round cap.
    size_t changepoint_capped = 0;
    /// Sites demoted from reliable to unreliable for low phase coherence, and how many had passed
    /// the `reliability` gate -- the gap between the two criteria, counted directly.
    size_t demoted_incoherent = 0;
    /// Demotion rounds actually run, and chains that were still demoting when the cap hit -- a
    /// chain that never reaches a coherent fixed point is telling you something and must not be
    /// silently reported as converged.
    size_t coherence_rounds_run = 0;
    size_t coherence_unconverged = 0;
    /// Sites the backbone local search flipped away from the cascade's answer, and rounds it took.
    size_t backbone_flips = 0;
    size_t backbone_rounds = 0;
    /// Sites kept out of the backbone by triangle consistency, and the triangles scored.
    size_t triangle_excluded = 0;
    size_t triangles_scored = 0;
    size_t triangles_open = 0;
    /// Sites where the lookback vote disagreed with the adjacent link the plain cascade would have
    /// obeyed -- the bad links it caught, counted directly.
    size_t lookback_overruled = 0;
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
