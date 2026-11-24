#ifndef VG_READ_CONSENSUS_STEP_HPP_INCLUDED
#define VG_READ_CONSENSUS_STEP_HPP_INCLUDED

/** \file
 * The step of the read-likelihood caller that gives a snarl, or a run of neighbouring snarls of
 * one chain, an allele built from its reads.
 */

#include <cstddef>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "allele_splicer.hpp"
#include "genotype_linker.hpp"
#include "phase_table.hpp"
#include "read_consensus.hpp"
#include "read_likelihood_caller.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"
#include "tree_genotyper.hpp"

namespace vg::multipass {

using namespace std;

/// The settings of the read-consensus step.
struct ReadConsensusOptions {
    /// Whether the step runs.
    bool enabled = false;
    /// The fewest bases by which a read differs from the sequence called on its strand for the
    /// read to disagree with it: a difference in length of that many bases, or that many
    /// differing bases in an alignment of the two that are also `min_diff_share` of the longer.
    size_t min_diff = 15;
    /// The fewest reads of a strand that disagree with its called sequence for the strand to be
    /// a candidate. It is at most `min_spanning`.
    size_t min_reads = 3;
    /// A consensus becomes an allele where it differs from the called sequence by more bases
    /// than this.
    size_t min_change = 10;

    /// The share of the longer of a read and the called sequence that their differing bases
    /// are, at least, where a read of nearly the called length disagrees.
    static constexpr double min_diff_share = 0.08;
    /// The least share of a strand's spanning reads that disagree, for the strand to be a
    /// candidate.
    static constexpr double min_share = 0.25;
    /// The fewest reads of a strand that span the site, for the strand to be a candidate.
    static constexpr size_t min_spanning = 6;
    /// A run grows from a site where `min_reads` of its spanning reads have this many bases
    /// more between its bounds than its longest called sequence, or this many fewer than its
    /// shortest, as it does from a site with a candidate strand.
    static constexpr size_t run_site_diff = 3;
    /// The most bases that lie between two neighbouring sites of a run: their bound nodes, and
    /// any nodes between those. Two sites whose reads have a gap of `min_diff` bases inside
    /// the node of a bound that faces the other are neighbours across more.
    static constexpr size_t run_max_junction = 32;
    /// The most sites a run has.
    static constexpr size_t run_max_members = 8;
};

/// What the read-consensus step did, for its report.
struct ReadConsensusCounts {
    /// Snarls with a phase and a read-likelihood call, which the step looked at.
    size_t sites = 0;
    /// Those whose reads were fetched and compared with the called sequences, on their own.
    size_t tested = 0;
    /// Strands that were candidates at one of those.
    size_t candidate_strands = 0;
    /// Runs treated: their consensus was cut and their member sites genotyped again, or they
    /// became one site.
    size_t runs_treated = 0;
    /// Treated runs that became one site.
    size_t runs_merged = 0;
    /// Sites that gained an allele, on their own or as a member of a treated run, and sites
    /// made of a run.
    size_t sites_changed = 0;
    /// Those whose chosen genotype has a consensus allele, after the last linkage pass.
    size_t sites_called = 0;
    /// Snarls left out because a called sequence has more bases between the bounds than the
    /// scorer aligns, and tested snarls where an allele has, or where the reads and alleles need
    /// more alignments than the scorer makes.
    size_t over_cap = 0;
    /// Candidate strands left alone because a site nested in theirs gained an allele for the
    /// strand.
    size_t covered_by_nested = 0;
};

/**
 * The step that adds to a snarl, or to the snarls of a run, the consensus of the reads that
 * disagree with what is called there, and genotypes them again.
 *
 * It runs after the first linkage pass and read phasing, on the staged snarls that have a phase.
 * At each, a strand's reads are the reads that span the snarl from bound to bound and that the
 * read strand table puts on the strand. Where the table gives a strand to fewer than half of
 * the spanning reads, all of them are taken to be on strand 0. A snarl with one strand has all
 * of them on it.
 *
 * A read disagrees with a strand where its bases between the bounds differ from the strand's
 * called sequence there (`AlleleSplicer::strand_sequence`) by `min_diff` bases or more (see
 * `ReadConsensusOptions::min_diff`). A strand is a candidate where at least `min_spanning` of its
 * reads span the snarl, and at least `min_reads` of them, and at least `min_share` of them,
 * disagree. A candidate strand gets the consensus of its spanning reads, with its called
 * sequence as the backbone (`read_consensus`). A consensus more than `min_change` differing
 * bases from the called sequence becomes a sequence allele of the snarl
 * (`add_sequence_allele`), unless an allele of the snarl spells it already. Two strands with one
 * consensus add one allele.
 *
 * A strand is left alone where a snarl staged under its snarl gained an allele for the strand,
 * so that a difference is given to the innermost snarl that holds it.
 *
 * A snarl that gained an allele is genotyped again from its reads, which the calculator aligns
 * to every allele, and filed again in the linkage model in place of its entry. The caller then
 * runs the linkage pass and read phasing again.
 *
 * Before any read is fetched, a snarl is left out where no read can disagree: where fewer than
 * `min_spanning` reads span it, or where every called sequence has fewer than `min_diff` bases
 * between the bounds and fewer than `min_reads` spanning reads have that many. The direct pass
 * keeps the numbers this needs with each snarl's call (`AlleleLikelihoodParams::span_rank`).
 *
 * ## Runs
 *
 * A mapper places a difference in a repeat in one of several neighbouring snarls, and not in
 * the same one for every read, so no one snarl need hold `min_diff` bases of it. A run is two or
 * more snarls that follow each other in one chain, at any level: each shares its bound node
 * with the next, or is joined to it by nodes that every walk of the chain takes. The junction
 * of two neighbours is the bases of those nodes, bound nodes included. A run is read from its
 * first member's start bound to its last member's end bound, and the same test, consensus and
 * threshold are applied to it as to one snarl: a strand's called sequence across the run is its
 * called sequence at each member with the junctions between, and a read's bases are those it
 * has between the run's two outer bounds.
 *
 * A snarl can be in a run where the single-snarl test gives it no allele, it has a called
 * sequence of no more bases than the scorer aligns, and at least `min_spanning` reads span it.
 * Two such snarls are neighbours where they follow each other in a chain, one phase set orders
 * the strands of both, and their junction has at most `run_max_junction` bases or holds a gap:
 * `min_reads` of the reads that span one of the two have a gap of `min_diff` bases or more
 * inside the node of its bound that faces the other (`SpanSummary::start_gap`, `end_gap`). A
 * stretch is a longest sequence of neighbours.
 *
 * A snarl's signal is the number of bases by which `min_reads` of its spanning reads have more
 * between the bounds than its longest called sequence, or fewer than its shortest, whichever is
 * more, or the length of the gap its junction with a neighbour holds, if that is more; the
 * direct pass's span summary gives both. A seed is a snarl of a stretch with a
 * candidate strand, or with a signal of `run_site_diff` or more. Runs are grown in a stretch
 * that has a seed and either a snarl with a candidate strand or signals that sum to `min_diff`
 * or more, so that reads are fetched for few stretches.
 *
 * A run grows from each seed in turn, in the stretch's order, among the snarls of the stretch
 * that are in no treated run. It takes the neighbour with the greater signal, or on a tie the
 * one on the side where it has fewer members, the later side if those tie too, and is tested:
 * its reads are fetched and compared as above. It stops at the first run that is treated, so
 * the run is the smallest around the seed that is. It also stops where the run can neither be
 * cut nor become one site, is over a limit of the scorer, has `min_spanning` reads on no strand, has `run_max_members`
 * members, or has no neighbour that keeps the bases between its bounds within what the scorer
 * aligns. A treated run then takes in one more neighbour, the later one first, for as long as
 * that leaves it treated and gives its treated strands more reads that disagree with what is
 * called: those are reads whose difference the mapper put in that neighbour. So a snarl is in
 * at most one treated run.
 *
 * The consensus of a candidate strand across a run is cut into the members: it is aligned to
 * the strand's called sequence, and each junction must be matched whole, base for base. Where
 * the best alignment does not match one, the two are aligned again with the junctions pinned
 * (`align_global_pinned`), and that alignment stands if it costs no more than the best plus the
 * cost of opening one long gap for each junction the best did not match. A run with no such
 * alignment becomes one site (below). Each member of a run that was cut has the consensus's
 * bases between its bounds, which become a sequence allele of it unless one of its alleles
 * spells them.
 *
 * Every read that spans a run that was cut is then placed in the members in one way. It is aligned to each
 * sequence known across the run, the strands' called sequences and the consensuses, and cut at
 * the junctions of the one it fits best. Each member of the run is genotyped again with those
 * bases for those reads (`ReadInteriors`), whether or not it gained an allele, and filed again.
 * A read that does not span the run is scored as at any other snarl.
 *
 * A run is left alone on a strand where a snarl staged under one of its members gained an
 * allele for the strand, and a member of a treated run counts as having gained one for the
 * strands the run was treated on.
 *
 * ## A run as one site
 *
 * A member cannot change the bases of a junction, so a run with a consensus that cannot be cut
 * becomes one site in place of its members: a merged site of the overlay of sites
 * (`EditSnarlDecompositionOverlay::add_merged_site`), from the first member's start bound to
 * the last member's end bound. `TreeGenotyper::merged_site` makes and genotypes it. Its
 * alleles are sequences across the run: the walk that takes each member's reference allele,
 * which is its reference allele; what each strand is called as across the run, as the walk
 * that takes the strand's allele at each member where that walk spells it, and as a sequence
 * allele otherwise; and each strand's consensus, whole, as a sequence allele. Alleles that
 * spell the same bases are one.
 *
 * The site is staged in the place of the run's first member, and filed in the linkage model in
 * place of that member's entry. The other members, every site staged under a member, and the
 * edit sites inside the nodes of the junctions are removed from the staged sites and from the
 * linkage model (`StagedSiteTable::remove_with_nested`), since the site's alleles spell what
 * was called at them.
 *
 * A run does not become one site, and is left alone on the strands whose consensus cannot be
 * cut, where its members are not placed alike in the nesting tree and on the reference, where
 * a site staged under a member gained an allele, or where a change is called at an edit site
 * inside a node of a junction, which the strands' called sequences across the run do not have.
 *
 * The reads are fetched with `SiteReadSource::for_each_read_uncounted`, both to compare them and
 * to genotype again, so the step adds nothing to the source's tallies of read starts or to its
 * off-panel candidates.
 */
class ReadConsensusStep {
public:
    /// The objects the step reads and writes. None of them is owned.
    struct Parts {
        const PathPositionHandleGraph* graph = nullptr;
        /// The overlay of sites that the staged sites' net handles are in.
        const EditSnarlDecompositionOverlay* sites = nullptr;
        /// The genotyper that called the sites, whose calculator fetches their reads.
        ReadLikelihoodSnarlCaller* genotyper = nullptr;
        /// Genotypes a site again with its added alleles.
        const TreeGenotyper* tree_genotyper = nullptr;
        /// Files a site in the linkage model again.
        const GenotypeLinker* linker = nullptr;
        StagedSiteTable* staged_sites = nullptr;
        /// The phase of each site. The step fills and then empties its render phases.
        PhaseTable* phase_table = nullptr;
        /// The strand of each read, with its sums built (`ReadStrandTable::build_lambda`).
        const ReadStrandTable* read_strands = nullptr;
    };

    /// Set the step's settings.
    void set_options(const ReadConsensusOptions& new_options) { options = new_options; }

    /// Get the step's settings.
    const ReadConsensusOptions& get_options() const { return options; }

    /// Return whether the step runs.
    bool enabled() const { return options.enabled; }

    /**
     * Run the step over the staged sites.
     *
     * Never call this inside a parallel region.
     *
     * @return The number of sites that gained an allele and were genotyped again.
     */
    size_t run(const Parts& parts);

    /// Count the sites that gained an allele, and the sites made of a run, whose chosen genotype
    /// has a consensus allele, once the last linkage pass has chosen the genotypes.
    void count_called(const Parts& parts);

    /// Write the step's one-line report.
    void report(ostream& out) const;

    /// Get what the step did.
    const ReadConsensusCounts& counts() const { return counted; }

private:
    /// What the step knows of one site on its own.
    struct Alone;
    /// A run of neighbouring sites of one chain, or one site.
    struct Run;
    /// What the step found at one run.
    struct RunFinding;

    /// Grow and test the runs of two or more among `snarls`, of each of which `alone` says what
    /// is known on its own. Add to `runs` each run that is treated, and to `run_findings` what
    /// was found at it.
    void grow_runs(const Parts& parts, const AlleleSplicer& splicer,
                   const vector<StagedSite*>& snarls, const vector<Alone>& alone,
                   vector<Run>& runs, vector<RunFinding>& run_findings);

    /// Compare the reads that span `run` with each strand's called sequence across it, build
    /// and cut the consensus of each candidate strand, and, for a run of two or more, find the
    /// bases every spanning read has between the bounds of each member site.
    RunFinding examine_run(const Parts& parts, const AlleleSplicer& splicer,
                           const vector<StagedSite*>& snarls, const Run& run) const;

    /// Return whether `run` can become one site, and note in `found` the edit sites between its
    /// members. It cannot where a change is called at one of them.
    bool can_merge(const Parts& parts, const vector<StagedSite*>& snarls, const Run& run,
                   RunFinding& found) const;

    /// Make the one site that `run` becomes, with what each strand is called as across the run
    /// and each consensus in `found` as its alleles, and genotype it.
    ///
    /// @return Nothing if no consensus is a new allele, or the site could not be genotyped.
    optional<StagedSite> merge_run(const Parts& parts, const AlleleSplicer& splicer,
                                   const vector<StagedSite*>& snarls, const Run& run,
                                   const RunFinding& found) const;

    /// Return what the span summary of `site` says of it: whether a read of it can disagree
    /// with a called sequence, and what it adds to the signal of a run.
    Alone measure(const Parts& parts, const AlleleSplicer& splicer, const StagedSite& site) const;

    /// Return whether a read with the bases `read` between the bounds disagrees with the called
    /// sequence `called`.
    bool disagrees(const string& read, const string& called) const;

    ReadConsensusOptions options;
    ReadConsensusCounts counted;
    /// The record keys of the sites that gained an allele, and of the sites made of a run.
    vector<size_t> changed;
    /// The staged edit sites inside each node, by node ID, while the step runs.
    unordered_map<nid_t, vector<const StagedSite*>> edit_sites_on;
    /// The most bases between the bounds, and the most alignments at one snarl, that the
    /// calculator takes at a snarl with a sequence allele
    /// (`AlleleLikelihoodCalculator::sequence_allele_limits`), while the step runs.
    size_t max_bases = 0;
    size_t max_alignments = 0;
};

}

#endif
