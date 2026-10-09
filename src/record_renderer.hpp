#ifndef VG_RECORD_RENDERER_HPP_INCLUDED
#define VG_RECORD_RENDERER_HPP_INCLUDED

#include <functional>
#include <vector>

#include "anchor.hpp"
#include "genotype_linker.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg {
namespace multipass {

using namespace std;

/**
 * The render: builds each staged site's records once, from its chosen genotype, after every
 * linkage pass and phasing round is done, and collects each site's anchors as it goes.
 */
class RecordRenderer {
public:
    /// Writes one staged site's records at `genotype`, through the caller that staged it.
    using LineWriter = function<void(const StagedSite& site, const vector<int>& genotype)>;

    /// Where to read a site's name.
    void configure(SiteReader reader);

    /// Hand every nested site that gets a line of its own to `staged`'s render queues, then write
    /// each queued site's records with `write_line`, a queue per thread, at the genotype `linker`
    /// chose. Each site's anchors go to `anchors`, unless it is null: just before its line, or, for
    /// a nested site with no line, before the hand-off. `phases` must be frozen for the render,
    /// and `strands` must hold each read's strand log-odds. Reports what it held back and how many
    /// records it rendered when `show_progress` is set.
    void render(StagedSiteTable& staged, const PhaseTable& phases, const ReadStrandTable& strands,
                const GenotypeLinker& linker, AnchorCollector* anchors,
                const LineWriter& write_line, bool show_progress) const;

private:
    /// Turn a staged site into anchors, with the phase order, the haploid slot, the leaf test and
    /// the gqn derived from it. The genotype is a parameter because the render passes the chosen
    /// pair. A site with no reference position still gets anchors, since a pin is placed by node
    /// ID.
    void collect_anchors(const StagedSite& site, const vector<int>& genotype,
                         const PhaseTable& phases, const ReadStrandTable& strands,
                         const LinkageCollector* model, AnchorCollector& anchors) const;

    /// The gqn column's value for a site: the direct pass's `gq_fraction`, unless the linkage model
    /// changed the call, in which case the signed value recomputed for the chosen genotype. NaN,
    /// written as `.`, where there is no value: no gap to normalise, or a moved call whose margin
    /// cannot be recomputed.
    static double anchor_gqn(const StagedSite& site, const vector<int>& chosen,
                             const LinkageCollector* model);

    /// Collect anchors for the nested sites that get no line, then move the others to the render
    /// queues (see `StagedSiteTable::hand_off`).
    void hand_off(StagedSiteTable& staged, const PhaseTable& phases,
                  const ReadStrandTable& strands, const GenotypeLinker& linker,
                  AnchorCollector* anchors, bool show_progress) const;

    SiteReader reader;
};

}
}

#endif
