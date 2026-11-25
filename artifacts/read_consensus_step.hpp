#ifndef VG_READ_CONSENSUS_STEP_HPP_INCLUDED
#define VG_READ_CONSENSUS_STEP_HPP_INCLUDED

/** \file
 * The step of the read-likelihood caller that gives a snarl, or a run of neighbouring snarls of
 * one chain, an allele built from its reads.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
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
    /// The fewest nats that the reads of a strand gain by a consensus, in all, for it to
    /// become an allele. A read gains what it scores more against the consensus than against
    /// what the strand is called as, and no more than one read's say
    /// (`AlleleLikelihoodCalculator::read_clamp`), and loses likewise.
    double min_gain = 9.0;
    /// How a strand's consensus is built, and how it is aligned to the walk the strand is
    /// called as to become an allele.
    ReadConsensusParams consensus;

    /// Return the fewest reads that can gain `min_gain` nats, where one read's say is `clamp`.
    size_t reads_needed(double clamp) const {
        return clamp > 0.0 ? (size_t)std::ceil(min_gain / clamp)
                           : std::numeric_limits<size_t>::max();
    }

    /// The most bases that lie between two neighbouring sites of a run, in their bound nodes
    /// and any nodes between those, and the most that lie between two seeds of a run that
    /// follow each other, in those nodes and in the sites between that are no seeds.
    static constexpr size_t run_max_junction = 1000;
};

/// What the read-consensus step did, for its report.
struct ReadConsensusCounts {
    /// Snarls with a phase and a read-likelihood call, which the step looked at.
    size_t sites = 0;
    /// Those whose reads were fetched and compared with the called sequences, on their own.
    size_t tested = 0;
    /// Strands at one of those whose consensus differs from what the strand is called as.
    size_t candidate_strands = 0;
    /// Runs of several snarls that were treated: each is to become one site.
    size_t runs_treated = 0;
    /// Those that became one site.
    size_t runs_merged = 0;
    /// Snarls that gained an allele, and sites made of a run.
    size_t sites_changed = 0;
    /// Those whose chosen genotype has a consensus allele, after the last linkage pass.
    size_t sites_called = 0;
    /// Snarls left out because a called sequence has more bases between the bounds than the
    /// scorer aligns of a read, and tested snarls where a consensus has, or where the reads and
    /// alleles need more alignments than the scorer makes.
    size_t over_cap = 0;
    /// Consensus alleles of those sites entered as the walk their strand takes across the site,
    /// with changes to its steps.
    size_t threaded = 0;
    /// Those entered as their bases between the bounds, because the strand has no walk across
    /// the site.
    size_t unspelled = 0;
};

/**
 * The step that adds to a snarl, or to the snarls of a run, the consensus of the reads of a
 * strand where the reads gain by it, and genotypes them again.
 *
 * It runs after the first linkage pass and read phasing, on the staged snarls that have a phase.
 * At each, a strand's reads are the reads that span the snarl from bound to bound and that the
 * read strand table puts on the strand. Where the table gives a strand to fewer than half of
 * the spanning reads, all of them are taken to be on strand 0. A snarl with one strand has all
 * of them on it.
 *
 * One threshold decides what is built and what is kept: `min_gain`, in nats. One read's say is
 * the most nats by which the genotyper lets one read favour an allele
 * (`AlleleLikelihoodCalculator::read_clamp`). A read's gain by a sequence is what it scores
 * more against that sequence than against its strand's called sequence
 * (`AlleleSplicer::strand_sequence`), each scored as a read is against a sequence allele
 * (`AlleleLikelihoodCalculator::sequence_nats`). It is no more than one read's say, and no less
 * than that below zero. So `reads_needed` reads at least gain `min_gain` nats.
 *
 * A strand with `reads_needed` spanning reads gets their consensus, with its called sequence as
 * the backbone (`read_consensus`). The strand is a candidate where the consensus differs from
 * the called sequence. The consensus becomes an allele where the strand's reads gain
 * `min_gain` nats by it in all, unless an allele of the snarl spells it already. The snarl is
 * then a run of one, and becomes one site as a run of several does (below). Two strands with
 * one consensus add one allele.
 *
 * The caller then runs the linkage pass and read phasing again.
 *
 * Before any read is fetched, a snarl is left out where its reads cannot gain `min_gain` nats
 * by an allele. No allele explains a spanning read that, between the bounds, scores `min_gain`
 * nats or more below a sequence that spells it against the allele it fits best. The direct
 * pass counts those reads with each snarl's call (`SpanSummary::unexplained`). Each can gain
 * one read's say, and a snarl is tested where that sums to `min_gain`.
 *
 * ## Runs
 *
 * A mapper places a difference in a repeat in one of several neighbouring snarls, and not in
 * the same one for every read, so no one snarl need hold the difference in enough reads, and
 * the bases that two reads have between the bounds of one snarl cannot be compared. So the step
 * tests a run of snarls as it tests one snarl. A run is read from its first member's start
 * bound to its last member's end bound: a strand's called sequence across the run is its
 * called sequence at each member with the junctions between, and a read's bases are those it
 * has between the run's two outer bounds.
 *
 * Two snarls are neighbours where one follows the other in a chain, at any level: they share a
 * bound node, or only nodes that every walk of the chain takes lie between them. One phase set
 * orders the strands of both, `reads_needed` reads span each, and each has a called sequence of
 * no more bases than the scorer aligns of a read. Their junction is the bases of the nodes from
 * the one's end bound to the other's start bound, both included, which are at most
 * `run_max_junction`.
 *
 * A snarl's signal is one read's say for each spanning read that no allele explains. A read's
 * gap is a run of bases it inserts in a node or beside it, or a run of a node's bases it
 * lacks. Where `reads_needed` of the reads that span a snarl have a gap inside the node of a
 * bound that faces a neighbour, or in the nodes of their junction beyond it
 * (`SpanSummary::start_gap`, `end_gap`), and a gap of that length costs a read `min_gain` nats,
 * the signal of both snarls is the say of that many reads at least. A seed is a snarl with a
 * signal of `min_gain` or more.
 *
 * A run begins and ends at a seed, and holds every snarl between two of its seeds. Two seeds
 * that follow each other in a sequence of neighbours are in one run where at most
 * `run_max_junction` bases lie between them, in the nodes of the junctions and in the snarls
 * that are no seeds, and the run has no more bases between its outer bounds than the scorer
 * aligns of a read. The bases of a snarl are those of its longest called sequence. Where one
 * more seed would pass either number, the run ends and the next begins at that seed. A run of
 * two or more is tested, and is treated where a strand has a consensus that is to become an
 * allele and the run can become one site. A snarl in no treated run is treated as a run of
 * one where it has a consensus, so a snarl is in at most one treated run.
 *
 * ## A run as one site
 *
 * A treated run becomes one site in place of its members: a merged site of the overlay of sites
 * (`EditSnarlDecompositionOverlay::add_merged_site`), from the first member's start bound to
 * the last member's end bound. `TreeGenotyper::merged_site` makes and genotypes it. Its
 * alleles are sequences across the run: the walk that takes each member's reference allele,
 * which is its reference allele; what each strand is called as across the run; and each
 * strand's consensus. Each of the last two is the walk that takes the strand's allele at each
 * member, with changes to the bases of its steps where it spells other bases
 * (`add_walk_allele`). Alleles that spell the same bases are one. It is genotyped from its
 * reads, which the calculator aligns to every allele.
 *
 * The site is staged in the place of the run's first member, and filed in the linkage model in
 * place of that member's entry. A run of one keeps its member's ID. The other members, every
 * site staged under a member, and the
 * edit sites inside the nodes of the junctions are removed from the staged sites and from the
 * linkage model (`StagedSiteTable::remove_with_nested`), since the site's alleles spell what
 * was called at them.
 *
 * A run cannot become one site, and is left alone, where its members are not placed alike in
 * the nesting tree and on the reference, or where a change is called at an edit site inside a
 * node of a junction, which the strands' called sequences across the run do not have.
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

    /// Find and test the runs of two or more among `snarls`, of each of which `alone` says what
    /// is known on its own. Add to `runs` each run that is treated, to `run_findings` what was
    /// found at it, and to `merged` and `entered` what `treat_run` made of it.
    void find_runs(const Parts& parts, const AlleleSplicer& splicer,
                   const vector<StagedSite*>& snarls, const vector<Alone>& alone,
                   vector<Run>& runs, vector<RunFinding>& run_findings,
                   vector<optional<StagedSite>>& merged,
                   vector<ReadConsensusCounts>& entered) const;

    /// Test `run` (`examine_run`) and, where it is treated, make its one site (`merge_run`)
    /// into `merged`, from the reads fetched once, which this thread holds until it is done.
    RunFinding treat_run(const Parts& parts, const AlleleSplicer& splicer,
                         const vector<StagedSite*>& snarls, const Run& run,
                         optional<StagedSite>& merged, ReadConsensusCounts& entered) const;

    /// Fetch the reads that span `run`, build the consensus of each strand's, and score the
    /// reads against it and against the strand's called sequence across the run. The run is
    /// treated for a consensus that is to become an allele: of the snarl, for a run of one, and
    /// of the one site that a run of several becomes.
    RunFinding examine_run(const Parts& parts, const AlleleSplicer& splicer,
                           const vector<StagedSite*>& snarls, const Run& run) const;

    /// Return whether `run` can become one site, and note in `found` the edit sites between its
    /// members. It cannot where a change is called at one of them.
    bool can_merge(const Parts& parts, const vector<StagedSite*>& snarls, const Run& run,
                   RunFinding& found) const;

    /// Make the one site that `run` becomes, with what each strand is called as across the run
    /// and each consensus in `found` as its alleles, and genotype it. Count in `entered` how
    /// each consensus is entered.
    ///
    /// @return Nothing if no consensus is a new allele, or the site could not be genotyped.
    optional<StagedSite> merge_run(const Parts& parts, const AlleleSplicer& splicer,
                                   const vector<StagedSite*>& snarls, const Run& run,
                                   const RunFinding& found, ReadConsensusCounts& entered) const;

    /// Return what the span summary of `site` says of it: whether its reads can gain
    /// `min_gain` nats by an allele, and its signal in a run.
    Alone measure(const Parts& parts, const AlleleSplicer& splicer, const StagedSite& site) const;

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
    /// The most nats by which one read can favour one allele over another
    /// (`AlleleLikelihoodCalculator::read_clamp`), while the step runs.
    double clamp = 0.0;
};

}

#endif
