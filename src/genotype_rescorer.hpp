#ifndef VG_GENOTYPE_RESCORER_HPP_INCLUDED
#define VG_GENOTYPE_RESCORER_HPP_INCLUDED

#include <functional>
#include <string>

#include "regenotype.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg {
namespace multipass {

using namespace std;

/**
 * Re-genotyping from the phase (see regenotype.hpp): rewrites each staged site's genotype
 * likelihoods from the strand log-odds of its reads, which read phasing gives them, so that the
 * linkage pass can choose the genotypes again from the result.
 */
class GenotypeRescorer {
public:
    /// Turn re-genotyping on or off, with `params`. `passes` is the most linkage passes; with 1
    /// the correction is computed and reported but not applied. `ledger`, if not empty, is where
    /// to write one line per site whose best genotype the correction changes.
    void configure(bool on, const RegenotypeParams& params, size_t passes, const string& ledger);

    /// Where to read what a staged site does not hold: the genotyper, which recomputes GQ, and a
    /// site's ID, for the ledger. Needed before `rescore`.
    void set_site_reader(SiteReader reader);

    /// Whether re-genotyping is on.
    bool enabled() const { return on; }
    /// See `configure`.
    size_t passes() const { return max_passes; }
    const RegenotypeParams& params() const { return regenotype_params; }

    /// The counters of the last `rescore`, which hold the temper and ceiling it used.
    const RegenotypeCounters& last_counters() const { return counters; }

    /// Re-score every staged site's genotype likelihoods with the reads' phase, and report the
    /// correction. `phases` and `strands` are as read phasing last left them, and `phase_set_id`
    /// numbers a contig's phase set for the run. The temper is fitted on the first round and kept
    /// in `fit` for the later ones.
    ///
    /// With more than one pass, the corrected likelihoods replace each site's own, always
    /// corrected from the direct pass's, and GQ is recomputed from them where the best genotype
    /// changed and is the direct pass's elsewhere. The answer at the other ploidy is corrected the
    /// same way. With one pass, nothing is kept.
    ///
    /// Returns true if any site's corrected best genotype differs from its called one. Does
    /// nothing, and returns false, when re-genotyping is off or `phases` is empty. The calibration
    /// table is reported only with `show_calibration`.
    bool rescore(StagedSiteTable& staged, const PhaseTable& phases, const ReadStrandTable& strands,
                 TemperFit& fit,
                 const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
                 bool show_calibration);

private:
    bool on = false;
    RegenotypeParams regenotype_params;
    size_t max_passes = 2;
    string ledger_path;
    SiteReader reader;
    RegenotypeCounters counters;
};

}
}

#endif
