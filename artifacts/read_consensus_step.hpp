#ifndef VG_READ_CONSENSUS_STEP_HPP_INCLUDED
#define VG_READ_CONSENSUS_STEP_HPP_INCLUDED

/** \file
 * The step of the read-likelihood caller that gives a snarl an allele built from its reads.
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
    /// The share of the longer of a read and the called sequence that their differing bases
    /// are, at least, where a read of nearly the called length disagrees.
    double min_diff_share = 0.08;
    /// The fewest reads of a strand that disagree with its called sequence for the strand to be
    /// a candidate.
    size_t min_reads = 3;
    /// The least share of a strand's spanning reads that disagree, for the strand to be a
    /// candidate.
    double min_share = 0.25;
    /// The fewest reads of a strand that span the site, for the strand to be a candidate.
    size_t min_spanning = 6;
    /// A consensus becomes an allele where it differs from the called sequence by more bases
    /// than this.
    size_t min_change = 10;
    /// The most reads a consensus is built from.
    size_t max_reads = 20;
    /// Whether a site is tested only where its span summary says a read can disagree. Without
    /// it every site is tested, which measures what the summary leaves out.
    bool prefilter = true;
};

/// What the read-consensus step did, for its report.
struct ReadConsensusCounts {
    /// Snarls with a phase and a read-likelihood call, which the step looked at.
    size_t sites = 0;
    /// Those whose span summary says a read can disagree with a called sequence.
    size_t prefiltered = 0;
    /// Those whose reads were fetched and compared with the called sequences.
    size_t tested = 0;
    /// Strands that were candidates.
    size_t candidate_strands = 0;
    /// The sites those strands are at.
    size_t candidate_sites = 0;
    /// The candidate sites the span summary would have left out. It is 0 with the prefilter
    /// on, which tests no such site.
    size_t candidates_outside_prefilter = 0;
    /// Consensuses built.
    size_t built = 0;
    /// Alleles added.
    size_t added = 0;
    /// Sites that gained an allele.
    size_t sites_with_added = 0;
    /// Sites whose chosen genotype has an added allele, after the last linkage pass.
    size_t sites_called = 0;
    /// Sites not tested because a called sequence or an allele has more bases between the
    /// bounds than the scorer aligns.
    size_t over_base_cap = 0;
    /// Candidate sites left alone because their reads and alleles need more alignments than
    /// the scorer makes.
    size_t over_alignment_cap = 0;
    /// Candidate strands left alone because a site nested in theirs gained an allele for the
    /// strand.
    size_t covered_by_nested = 0;
    /// Candidate strands with no consensus, for want of reads or of an alignment within the
    /// aligner's limits.
    size_t no_consensus = 0;
    /// Consensuses within `min_change` bases of the called sequence.
    size_t unchanged = 0;
    /// Consensuses that an allele of the site already spells.
    size_t already_allele = 0;
    /// Sites that could not be genotyped with the added alleles, which keep their call.
    size_t regenotype_failed = 0;
};

/**
 * The step that adds to a snarl the consensus of the reads that disagree with what is called
 * there, and genotypes the snarl again.
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
 * keeps the two numbers this needs with each snarl's call (`AlleleLikelihoodParams::span_rank`).
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
        /// The most bases between the bounds, and the most alignments at one snarl, that the
        /// scorer takes at a snarl with a sequence allele.
        size_t max_bases = 10000;
        size_t max_alignments = 2000;
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

    /// Return whether allele `allele` of the staged site `site` is one the step adds: an allele
    /// of a snarl that its walk does not spell.
    static bool is_consensus_allele(const EditSnarlDecompositionOverlay& sites,
                                    const StagedSite& site, int allele);

    /// Count the sites that gained an allele whose chosen genotype has one, once the last
    /// linkage pass has chosen the genotypes.
    void count_called(const Parts& parts);

    /// Write the step's one-line report.
    void report(ostream& out) const;

    /// Get what the step did.
    const ReadConsensusCounts& counts() const { return counted; }

private:
    /// What the step found at one site.
    struct Finding;

    /// Compare the reads of `site` with its called sequences, and build the consensus of each
    /// candidate strand.
    Finding examine(const Parts& parts, const AlleleSplicer& splicer,
                    const StagedSite& site) const;

    /// Return whether a read with the bases `read` between the bounds disagrees with the called
    /// sequence `called`.
    bool disagrees(const string& read, const string& called) const;

    ReadConsensusOptions options;
    ReadConsensusCounts counted;
    /// The record keys of the sites that gained an allele.
    vector<size_t> changed;
};

}

#endif
