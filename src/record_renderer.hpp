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
 * A renderer that writes the record of each staged site once, from its chosen genotype.
 *
 * It runs after every linkage pass and phasing round is done, and collects the anchors of each
 * site as it goes.
 */
class RecordRenderer {
public:
    /// Writes one staged site's records at `genotype`, through the caller that staged it.
    using LineWriter = function<void(const StagedSite& site, const vector<int>& genotype)>;

    /// Set where to read the name of a site.
    void configure(SiteReader site_reader);

    /**
     * Write the records of every staged site at the genotype `linker` chose.
     *
     * Every nested site that gets a line of its own goes to the render queues of `staged`, and
     * `write_line` writes the queued sites, one queue per thread. `before_lines`, if set, runs
     * once before the first line.
     *
     * The anchors of each site go to `anchors`, unless it is null: just before its line, or for a
     * nested site with no line, before the hand-off. `phases` must be frozen for the render, and
     * `strands` must hold the strand log-odds of each read.
     *
     * With `show_progress`, report what was held back and how many records were written.
     */
    void render(StagedSiteTable& staged, const PhaseTable& phases, const ReadStrandTable& strands,
                const GenotypeLinker& linker, AnchorCollector* anchors,
                const LineWriter& write_line, bool show_progress,
                const function<void()>& before_lines = nullptr) const;

private:
    /// Collect the anchors of a staged site at `genotype`.
    ///
    /// A site with no reference position still gets anchors, since an anchor is placed by node ID.
    /// A site with no read evidence of its own gets none.
    void collect_anchors(const StagedSite& site, const vector<int>& genotype,
                         const PhaseTable& phases, const ReadStrandTable& strands,
                         const LinkageCollector* model, AnchorCollector& anchors) const;

    /// Return the value of the gqn column of the anchors of a site.
    ///
    /// Where the linkage model left the call alone, this is `gq_fraction` from the direct pass. Where
    /// the model changed the call, `gq_fraction` describes the abandoned genotype, so this is the
    /// signed margin of the chosen genotype, scaled as GQN in the VCF is (see
    /// `ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype`). It is NaN, written as `.`,
    /// where there was no gap to normalise or the margin cannot be recomputed.
    static double anchor_gqn(const StagedSite& site, const vector<int>& chosen,
                             const LinkageCollector* model);

    /// Collect anchors for the nested sites that get no line, then move the others to the render
    /// queues (see `StagedSiteTable::hand_off`).
    ///
    /// A dropped site, whose sample has no copy of it, gets neither.
    void hand_off(StagedSiteTable& staged, const PhaseTable& phases,
                  const ReadStrandTable& strands, const GenotypeLinker& linker,
                  AnchorCollector* anchors, bool show_progress) const;

    /// Where to read the name of a site.
    /// TODO: Nothing reads this. Drop it and `configure`, with the call in multipass_caller.cpp?
    SiteReader reader;
};

}

#endif
