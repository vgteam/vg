#ifndef VG_READ_PHASER_HPP_INCLUDED
#define VG_READ_PHASER_HPP_INCLUDED

#include <functional>
#include <string>

#include "read_phasing.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg {

using namespace std;

/**
 * Read phasing (--read-phasing; see read_phasing.hpp): decides which strand carries which allele
 * at each diploid heterozygous site from the reads that span several sites, and swaps the strands
 * of the linkage model's phase where the reads disagree with it. Genotypes are not changed.
 */
class ReadPhaser {
public:
    /// Turn read phasing on or off, with `params`.
    void configure(bool on, const ReadPhasingParams& params);

    /// Whether read phasing is on.
    bool enabled() const { return on; }

    /// Phase `phases`, which a linkage pass has just filled, from the reads of `staged`, and report
    /// what it did. Each diploid heterozygous site's per-read evidence, reduced to its phased pair,
    /// goes to `strands.sites()`, and the sites whose pair the reads reverse to `strands.flips()`;
    /// the swaps are then made in `phases`, and carried down to the nested strands under them.
    /// `phase_set_id` numbers a contig's phase set for the run.
    ///
    /// Does nothing when read phasing is off or `phases` is empty, and leaves `strands.flips()` as
    /// they were when no site has evidence, so that a later round can keep an earlier one's.
    void phase(StagedSiteTable& staged, PhaseTable& phases, ReadStrandTable& strands,
               const function<size_t(const string& contig, size_t phase_set)>& phase_set_id);

private:
    bool on = false;
    ReadPhasingParams params;
    /// What the last `phase` did, for its report.
    ReadPhasingCounters counters;
};

}

#endif
