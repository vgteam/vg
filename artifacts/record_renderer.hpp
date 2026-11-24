#ifndef VG_RECORD_RENDERER_HPP_INCLUDED
#define VG_RECORD_RENDERER_HPP_INCLUDED

#include <functional>
#include <vector>

#include "anchor.hpp"
#include "genotype_linker.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * The writer that turns staged sites into records, reading phases, strands and chosen genotypes
 * from the tables it was built with.
 *
 * It writes the record of each staged site once, from its chosen genotype. It runs after every
 * linkage pass and phasing round is done, and collects the anchors of each site as it goes, if
 * the anchor collector it was built with is enabled.
 */
class RecordRenderer {
public:
    /// Writes one staged site's records at `genotype`, through the caller that staged it. It
    /// may keep on the site what writing its records needs.
    using LineWriter = function<void(StagedSite& site, const vector<int>& genotype)>;

    /**
     * Make a renderer that reads from the given tables and collects anchors into `anchors`.
     *
     * Nothing is owned. The tables are read only during `render`, so they can still be filled
     * after construction.
     *
     * @param phases The phase of each site. It must be frozen for the render before `render`.
     * @param strands The strand log-odds of each read, for the anchors.
     * @param linker Gives the chosen genotype of each site.
     * @param anchors Collects anchors during `render` if it is enabled then.
     */
    RecordRenderer(const PhaseTable& phases, const ReadStrandTable& strands,
                   const GenotypeLinker& linker, AnchorCollector& anchors);

    /**
     * Write the records of every staged site at its chosen genotype.
     *
     * Every nested site that gets a line of its own goes to the render queues of `staged`, and
     * `write_line` writes the queued sites, one queue per thread. `before_lines`, if set, runs
     * once before the first line.
     *
     * Each site's anchors are collected just before its line, or, for a nested site with no line,
     * before the hand-off.
     *
     * With `show_progress`, report what was held back and how many records were written.
     */
    void render(StagedSiteTable& staged, const LineWriter& write_line, bool show_progress,
                const function<void()>& before_lines = nullptr) const;

private:
    /// Collect the anchors of a staged site at `genotype`.
    ///
    /// A site with no reference position still gets anchors, since an anchor is placed by node ID.
    /// A site with no read evidence of its own gets none.
    void collect_anchors(const StagedSite& site, const vector<int>& genotype) const;

    /// Return the value of the gqn column of the anchors of a site.
    ///
    /// Where the linkage model left the call alone, this is `gq_fraction` from the direct pass. Where
    /// the model changed the call, `gq_fraction` describes the abandoned genotype, so this is the
    /// signed phred margin of the chosen genotype over the best other genotype of the written
    /// alleles. It is divided by the achievable gap of the direct pass, multiplied by its explained
    /// share and clamped to [-1, 1], so that it agrees with the GQN that
    /// `ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype` writes. It is NaN, written
    /// as `.`, where there was no gap to normalise or the margin cannot be recomputed.
    static double anchor_gqn(const StagedSite& site, const vector<int>& chosen,
                             const LinkageCollector* model);

    /// Collect anchors for the nested sites that get no line, then move the others to the render
    /// queues (see `StagedSiteTable::hand_off`).
    ///
    /// A dropped site, whose sample has no copy of it, gets neither.
    void hand_off(StagedSiteTable& staged, bool show_progress) const;

    /// Gives the phase-ordered genotype and haploid slot of each site, for the anchors.
    const PhaseTable& phases;
    /// Gives each read's strand log-odds, for the anchors.
    const ReadStrandTable& strands;
    /// Gives the chosen genotype of each site, and the linkage model that chose it.
    const GenotypeLinker& linker;
    /// Collects the anchors of each site, when it is enabled.
    AnchorCollector& anchors;
};

}

#endif
