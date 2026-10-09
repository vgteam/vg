#ifndef VG_GENOTYPE_LINKER_HPP_INCLUDED
#define VG_GENOTYPE_LINKER_HPP_INCLUDED

#include <functional>
#include <string>
#include <vector>

#include "handle.hpp"
#include "linkage_model.hpp"
#include "snarl_caller.hpp"
#include "child_placer.hpp"
#include "panel_lookup.hpp"
#include "phase_table.hpp"
#include "staged_site.hpp"

namespace vg {
namespace multipass {

using namespace std;

/// A site's place on the reference: a contig, and an offset along it. Positions on one contig
/// order its sites, and the difference of two positions is the distance in bases between the
/// sites, except where a position is a stand-in (see `position`). The linkage model relies on
/// both: it builds its top-level linkage chains from the sites of one contig, orders a chain's
/// sites by position and measures the gaps between them from it, and names a phase set
/// (FORMAT/PS) by the position of the chain's first site.
struct SiteLocus {
    /// The contig as the VCF names it: the locus part of a PanSN path name, so `chr20` for
    /// `CHM13#0#chr20`, or the path name itself when it is not a PanSN name.
    string contig;
    /// A 0-based offset along the contig. For a site the reference path passes through, it is
    /// where the site's first boundary node starts on that path. It must not depend on which
    /// alleles the site's records carry, so it is not their POS, which trimming the alleles
    /// can move. For a site that no reference path passes through, it is a stand-in (see
    /// `GenotypeLinker::off_reference_site_locus`), and its difference from another position is a
    /// distance only for a site of the same chain, or for the chain's parent.
    size_t position = 0;
};

/**
 * The caller's side of the linkage model (`LinkageCollector`), which chooses genotypes jointly
 * along chains of sites from the haplotype panel.
 *
 * The direct pass files each genotyped site with `add`. The linkage pass, `link`, then has the
 * model choose the genotypes one level of the nesting tree at a time, parents before their
 * children, and revises each nested chain in the staged sites from its parent's chosen genotype:
 * its place along the parent, how many copies of it the sample has, and whether it exists at
 * all. Re-genotyping gives the model corrected likelihoods with `resync` and links again.
 *
 * Without a collector, nothing is filed or chosen, but `link` still revises the nested chains
 * from their parents' direct calls.
 */
class GenotypeLinker {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// Link through `collector`, which is not owned; null turns linkage off. `panel`, also not
    /// owned, gives each site's panel alleles.
    void configure(LinkageCollector* collector, const PanelLookup* panel);

    /// Where to read what a staged site does not hold: the graph, for a site's locus, its alleles'
    /// sequences and where its chains start, and the read-likelihood genotyper, whose GQ settings
    /// (the share and depth discounts) the quality inputs filed with each site follow. Needed
    /// before `add` or `link`.
    void set_site_reader(SiteReader reader);

    /// Whether there is a linkage model.
    bool enabled() const { return model != nullptr; }

    /// The linkage model, or null.
    LinkageCollector* collector() const { return model; }

    /// The locus of a site that the reference path `ref_path_name` passes through. `ref_offset`
    /// is added to every position along that path, to place the path on its contig.
    SiteLocus site_locus(const SiteBounds& site, const string& ref_path_name,
                         int ref_offset) const;

    /// The locus of a site that no reference path passes through, and that so has no position of
    /// its own. `ref_path_name` is the reference path through the site's nearest ancestor on a
    /// reference path, and `stand_in_position` is that ancestor's position plus how far along the
    /// ancestor's allele the site's chain starts (`StagedSite::position_from_parent`).
    static SiteLocus off_reference_site_locus(const string& ref_path_name,
                                              int64_t stand_in_position);

    /// File a site in the linkage model when the direct pass genotypes it, rather than when its
    /// line is written, since the linkage pass reads the model before any line is written. The
    /// written allele map and whether a line was written are supplied later, by
    /// `LinkageCollector::set_allele_map`. The direct pass does not call this for a retained chain
    /// with a reference path (see `NestingPlacement::retain_only`); the linkage pass files that
    /// chain if the sample carries it. Safe to call from several threads.
    ///
    /// `record_key` names the site (see `VCFOutputCaller::record_key_of`), and `placement` places
    /// it in the nesting tree. `ref_path_name` and `ref_offset` give the site's locus as for
    /// `site_locus`. `no_reference` marks a site that no reference path passes through:
    /// `ref_path_name` is then the reference path through the site's nearest ancestor on a
    /// reference path, `position_from_parent` is the site's stand-in position (see
    /// `off_reference_site_locus`), and `ref_offset` is not used. Otherwise
    /// `position_from_parent` is not used.
    ///
    /// Returns whether the site was filed: only a read-likelihood call (`score` not null) of one or
    /// two alleles, none of them missing, is. If it was and `panel_out` is given, the site's panel
    /// alleles are moved to `panel_out`, so that the staged site can keep them rather than look
    /// them up again.
    bool add(const SiteBounds& site, const vector<Traversal>& travs,
             const vector<int>& trav_genotype, const SiteScore* score,
             int ref_trav_idx, const string& ref_path_name, int ref_offset, size_t record_key,
             const NestingPlacement& placement, bool no_reference = false,
             int64_t position_from_parent = 0, vector<int>* panel_out = nullptr) const;

    /// The frequency exponent a site should decode with: `LinkageModel::Params::hp_prior` at a
    /// run-length site, or -1 for the model's own. Reads the traversals' sequences only when
    /// `hp_prior` is on.
    double freq_prior(const vector<Traversal>& travs, int ref_trav_idx) const;

    /// The genotype the linkage model chose for a staged site, or the direct pass's genotype
    /// where the model chose none.
    vector<int> chosen_genotype(const StagedSite& site) const;

    /// What one linkage pass did, for the report.
    struct PassCounts {
        /// The deepest level linked.
        size_t levels = 0;
        /// Nested chains filed again at a new ploidy, and filed for the first time.
        size_t revised = 0, gained = 0;
        /// Entries taken out of the model because the parent's chosen genotype does not carry
        /// their chain.
        size_t retracted = 0;
        /// Chains whose crossing mask the direct pass could not compute.
        size_t crossing_unknown = 0;
        /// Chains the model could not file again, for want of a compact allele space.
        size_t unspecifiable = 0;
        /// Chains dropped because no candidate traversal of the parent crosses them.
        size_t no_crossing = 0;
        /// Chains whose parent's chosen pair could not be read.
        size_t no_chosen = 0;
        /// Chains that cannot be built, so are left as they are.
        size_t unrenderable = 0;
        /// Chains left at a ploidy the direct pass never scored.
        size_t ploidy_unscored = 0;
    };

    /// The linkage pass: choose the genotypes one level at a time, and revise the nested chains of
    /// `sites` level by level. The model chooses each level's genotypes, adding their phase calls
    /// to `phases` when `keep_phase` is set; then each chain of the next level is placed along its
    /// parent's chosen allele and takes the ploidy the parent's chosen genotype gives it, from the
    /// answers the direct pass kept at both ploidies. A chain the parent does not carry is dropped
    /// with everything inside it. Each pass decides all of this again, so `phases` is cleared on a
    /// later pass. Does nothing to a site that is not nested. Never call it inside a parallel
    /// region.
    PassCounts link(StagedSiteTable& sites, PhaseTable& phases, bool keep_phase);

    /// Give the linkage model each staged site's likelihoods again, as re-genotyping corrected
    /// them, with the corrected best genotype as the called pair, and report how many it took.
    /// The genotypes are chosen again by the next `link`.
    void resync(StagedSiteTable& sites) const;

    /// Report the model's size, how many genotypes it moved and the time it took. Does nothing
    /// without a linkage model.
    void report() const;

private:
    /// Have the model choose the genotypes of one level, adding their phase calls to `calls`
    /// unless it is null. The phase calls accumulate across levels, since the model reads the
    /// earlier ones. `last` marks a linkage pass's final level, which builds the model's outputs
    /// from everything chosen.
    void resolve_level(size_t level, bool last, vector<PhaseCall>* calls);

    /// Drop the nested site `root` and its whole subtree: the chosen parent does not carry the
    /// chain, so the sample has no copy of it or of anything inside it. Returns how many entries
    /// were taken out of the model.
    size_t drop_subtree(StagedSiteTable& sites, size_t root);

    LinkageCollector* model = nullptr;
    const PanelLookup* lookup = nullptr;
    SiteReader reader;

    /// How many times the linkage pass has run.
    size_t passes_run = 0;
    /// Time spent in the linkage model, over every level of every pass, and how many genotypes
    /// the last pass moved off the direct calls, for the report.
    double total_seconds = 0.0;
    size_t total_moved = 0;
};

}
}

#endif
