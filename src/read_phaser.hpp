#ifndef VG_READ_PHASER_HPP_INCLUDED
#define VG_READ_PHASER_HPP_INCLUDED

#include <functional>
#include <string>

#include "read_phasing.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A phaser that orders the strands at heterozygous sites by the reads that span several sites.
 *
 * It decides which strand carries which allele at each diploid heterozygous site, and swaps the
 * strands of the linkage model's phase where the reads disagree with it (see read_phasing.hpp).
 * It does not change genotypes.
 */
class ReadPhaser {
public:
    /// Turn read phasing on or off, and set its parameters.
    void configure(bool enable, const ReadPhasingParams& settings);

    /// Return whether read phasing is on.
    bool enabled() const { return on; }

    /**
     * Re-phase `phases` from the reads of `staged`, and report what changed.
     *
     * The read evidence of each diploid heterozygous site goes to `strands.sites()`, and the
     * record keys of the sites whose strands the reads swap go to `strands.flips()`. The swaps
     * are then made in `phases`, and carried down to the nested sites.
     *
     * This does nothing when read phasing is off or `phases` is empty. When no site has read
     * evidence, it leaves `strands.flips()` as it was, so a later round keeps an earlier round's
     * swaps.
     *
     * @param phase_set_id Returns the run-wide ID of a phase set of a contig.
     */
    void phase(StagedSiteTable& staged, PhaseTable& phases, ReadStrandTable& strands,
               const function<size_t(const string& contig, size_t phase_set)>& phase_set_id);

private:
    /// Fill `sites` with the read evidence of each diploid heterozygous site, in record order.
    ///
    /// `phase_index` is `phases.index()`.
    void collect_sites(StagedSiteTable& staged, const PhaseTable& phases,
                       const unordered_map<size_t, size_t>& phase_index,
                       vector<PhaseSite>& sites,
                       const function<size_t(const string& contig, size_t phase_set)>&
                           phase_set_id) const;

    /// Write what the last `phase()` did to standard error.
    void report() const;

    /// Whether read phasing is on.
    bool on = false;
    /// The parameters of read phasing.
    ReadPhasingParams params;
    /// What the last `phase()` did.
    ReadPhasingCounters counters;
};

}

#endif
