#ifndef VG_GENOTYPE_LINKER_HPP_INCLUDED
#define VG_GENOTYPE_LINKER_HPP_INCLUDED

#include <functional>
#include <string>
#include <vector>

#include "../src/handle.hpp"
#include "linkage_model.hpp"
#include "../src/snarl_caller.hpp"
#include "descent_policy.hpp"
#include "panel_lookup.hpp"
#include "phase_table.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A linker that files the caller's sites in the linkage model and applies the model's choices.
 *
 * The linkage model (`LinkageCollector`) chooses genotypes jointly along runs of sites, from the
 * haplotype panel. The direct pass files each genotyped site with `file()`. The linkage pass,
 * `link()`, then has the model choose the genotypes one level of the snarl tree at a time,
 * parents before their children. After each level it revises each nested chain from the genotype
 * chosen for its parent: where the chain is along the parent, how many copies the sample has,
 * and whether it is there at all. Re-genotyping gives the model corrected likelihoods with
 * `resync()` and links again.
 *
 * Without a linkage model, nothing is filed or chosen, and each nested chain keeps the ploidy
 * the direct pass gave it.
 */
class GenotypeLinker {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// Set the linkage model, or null for none, and the panel lookup. Neither is owned.
    void configure(LinkageCollector* collector, const PanelLookup* panel);

    /// Set where to find the graph, the genotyper and the allele sequences. Call this before
    /// `file()` or `link()`.
    void set_site_reader(SiteReader site_reader);

    /// Return whether there is a linkage model.
    bool enabled() const { return model != nullptr; }

    /// Return the linkage model, or null.
    LinkageCollector* collector() const { return model; }

    /**
     * File a staged site in the linkage model, with its genotype as the called pair.
     *
     * The direct pass files each site when it is genotyped, because the linkage pass reads the
     * model before any record is written. `LinkageCollector::set_allele_map()` later supplies the
     * written alleles. The linkage pass files a nested chain again when it revises the chain's
     * ploidy.
     *
     * Only a read-likelihood call (`score` not null) of one or two alleles, none of them
     * missing, is filed. Filing fills the site's `panel_cache`. This is safe to call from
     * several threads on different sites. It adds an entry even where the key has a live one,
     * so the linkage pass must retract the old entry first.
     *
     * @param nested Whether the site's chain has one copy (see
     * `LinkageCollector::SiteContext::nested`). The direct pass gives
     * `NestingPlacement::one_copy`, and the linkage pass a ploidy of 1. These differ for a
     * haploid site at the top level.
     * @return Whether the site was filed.
     */
    bool file(StagedSite& site, bool nested) const;

    /// Return the frequency exponent a site decodes with, or -1 for the model's own.
    ///
    /// A run-length site takes `LinkageModel::Params::hp_prior` when that is on.
    double freq_prior(const CandidateAlleles& candidates) const;

    /// Return the genotype the linkage model chose for a staged site, or else its direct-pass
    /// genotype.
    vector<int> chosen_genotype(const StagedSite& site) const;

    /// What one linkage pass did, for the report.
    struct PassCounts {
        /// The deepest level linked.
        size_t levels = 0;
        /// Nested chains filed again at a new ploidy.
        size_t revised = 0;
        /// Nested chains filed for the first time.
        size_t gained = 0;
        /// Entries taken out of the model because the genotype chosen for the parent does not
        /// carry their chain.
        size_t retracted = 0;
        /// Chains whose crossing mask the direct pass could not compute.
        size_t crossing_unknown = 0;
        /// Chains the model could not file again, for want of a compact allele space.
        size_t unspecifiable = 0;
        /// Chains dropped because no candidate traversal of the parent crosses them.
        size_t no_crossing = 0;
        /// Chains for whose parent no chosen pair could be read.
        size_t no_chosen = 0;
        /// Chains left as they are because no record can be built from them.
        size_t unrenderable = 0;
        /// Chains left at a ploidy the direct pass never scored.
        size_t ploidy_unscored = 0;
    };

    /**
     * Run the linkage pass, choosing genotypes and revising nested chains one level at a time.
     *
     * The model chooses the genotypes of each level, adding their phase calls to `phases` if
     * `keep_phase` is set. Then each chain of the next level is placed along the allele chosen
     * for its parent, and takes the ploidy the parent's genotype gives it, from the answers the
     * direct pass kept at both ploidies. A chain the parent does not carry is dropped with
     * everything inside it. Each pass decides all of this again, so a later pass clears
     * `phases` first.
     *
     * Never call this inside a parallel region.
     */
    PassCounts link(StagedSiteTable& sites, PhaseTable& phases, bool keep_phase);

    /// Give the linkage model the likelihoods of each staged site as re-genotyping corrected them.
    ///
    /// The corrected best genotype becomes the called pair. The next `link()` chooses the
    /// genotypes again.
    void resync(StagedSiteTable& sites) const;

    /// Report the size of the model, how many genotypes it moved and the time it took.
    void report() const;

private:
    /// Have the model choose the genotypes of one level.
    ///
    /// The phase calls are added to `calls`, unless it is null. They accumulate across levels,
    /// since the model reads those of the parents. `last` marks the last level of a pass, which
    /// builds the outputs of the model from everything chosen.
    void resolve_level(size_t level, bool last, vector<PhaseCall>* calls);

    /// Index the phase calls of the parents of the nested sites at `level`, by record key.
    ///
    /// Where a parent has several, the last one wins.
    static unordered_map<size_t, const PhaseCall*> parent_phase_calls(
        const vector<StagedSite>& pending, const PhaseTable& phases, size_t level);

    /**
     * Revise one nested chain from the genotype chosen for its parent.
     *
     * The chain is placed along the parent and given the ploidy that the genotype of the parent
     * implies, or dropped with its subtree if the parent does not carry it.
     *
     * @param index The index of the chain in `sites.nested()`.
     * @param chosen The phase calls of the parents at this level.
     * @param crossing_finders The crossing finder of each parent traversal, filled as needed.
     */
    void revise_chain(StagedSiteTable& sites, size_t index,
                      const unordered_map<size_t, StagedSite*>& record_by_key,
                      const unordered_map<size_t, const PhaseCall*>& chosen,
                      unordered_map<const Traversal*, ChildCrossingFinder>& crossing_finders,
                      PassCounts& counts);

    /// Move a nested chain to where it starts along the allele chosen for its parent.
    ///
    /// @param crossing_finders The crossing finder of each parent traversal, by its address,
    /// filled as needed.
    void place_chain(StagedSite& pr, const StagedSite& parent,
                     unordered_map<const Traversal*, ChildCrossingFinder>& crossing_finders) const;

    /// Find the traversals chosen for the parent of a nested chain, the second -1 for a haploid
    /// parent.
    ///
    /// They come from the phase call of the parent if it has one, and else from its chosen
    /// genotype. They are traversals, which the crossing mask of the chain is indexed by.
    ///
    /// @return Whether a chosen traversal was found.
    bool chosen_pair_of_parent(const StagedSite& pr, const StagedSite* parent,
                               const unordered_map<size_t, const PhaseCall*>& chosen, int& first,
                               int& second) const;

    /// Return whether a record can be built from a staged site, with traversals and a genotype
    /// in range.
    static bool buildable(const StagedSite& pr);

    /// Give a nested chain its answer at ploidy `copies`, swapping in the answer the direct pass
    /// kept at the other ploidy if needed.
    ///
    /// @return Whether the chain has an answer at that ploidy.
    static bool adopt_ploidy(StagedSite& pr, int copies, PassCounts& counts);

    /// Compute the crossing masks of the children of a chain again, from its traversals.
    void recompute_child_crossings(StagedSiteTable& sites, const StagedSite& pr) const;

    /// Drop the nested site `root` and its whole subtree.
    ///
    /// The parent chosen does not carry the chain, so the sample has no copy of it or of
    /// anything inside it.
    ///
    /// @return How many entries were taken out of the model.
    size_t drop_subtree(StagedSiteTable& sites, size_t root);

    /// The linkage model, or null for none.
    LinkageCollector* model = nullptr;
    /// Looks up the panel alleles of each site.
    const PanelLookup* lookup = nullptr;
    /// Finds the graph, the genotyper and the allele sequences.
    SiteReader reader;

    /// How many times the linkage pass has run.
    size_t passes_run = 0;
    /// The seconds spent in the linkage model, over every level of every pass.
    double total_seconds = 0.0;
    /// How many genotypes the last pass moved off the direct calls.
    size_t total_moved = 0;
};

}

#endif
