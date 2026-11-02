#ifndef VG_GENOTYPE_RESCORER_HPP_INCLUDED
#define VG_GENOTYPE_RESCORER_HPP_INCLUDED

#include <functional>
#include <string>
#include <string_view>

#include "regenotype.hpp"
#include "phase_table.hpp"
#include "read_strand_table.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A rescorer that corrects each staged site's genotype likelihoods with the phase of its reads.
 *
 * Read phasing gives each read a strand log-odds. The rescorer rewrites each site's likelihoods
 * from them (see regenotype.hpp), so that the linkage pass can choose the genotypes again.
 */
class GenotypeRescorer {
public:
    /**
     * Turn re-genotyping on or off, and set its parameters.
     *
     * @param passes The most linkage passes. With 1, the correction is computed and reported
     * but not applied.
     * @param ledger The file to write one line to for each site whose best genotype the
     * correction changes, or empty for none.
     */
    void configure(bool enable, const RegenotypeParams& settings, size_t passes,
                   std::string_view ledger);

    /// Set where to find the genotyper, which recomputes GQ. Call this before `rescore()`.
    void set_site_reader(SiteReader site_reader);

    /// Return whether re-genotyping is on.
    bool enabled() const { return on; }
    /// Return the most linkage passes.
    size_t passes() const { return max_passes; }
    /// Return the parameters of re-genotyping.
    const RegenotypeParams& params() const { return regenotype_params; }

    /// Return the counters of the last `rescore()`, which hold the calibration it used.
    const RegenotypeCounters& last_counters() const { return counters; }

    /**
     * Correct every staged site's genotype likelihoods with the phase of its reads, and report
     * the correction.
     *
     * With more than one pass, the corrected likelihoods replace each site's own. Each round
     * corrects the direct pass's likelihoods, and recomputes GQ where the best genotype changed.
     * The answer at the other ploidy is corrected the same way. With one pass, the correction is
     * thrown away.
     *
     * Re-genotyping that is off, or an empty `phases`, makes no correction.
     *
     * @param phases The phase, as read phasing last left it.
     * @param strands The read evidence and strand swaps of read phasing.
     * @param fit The temper, its ceiling and its calibration table, fitted on the first round
     *     and reused by the later ones.
     * @param phase_set_id Returns the run-wide ID of a phase set of a contig.
     * @param show_calibration Whether to report the calibration table of the temper fit.
     * @return Whether any site's corrected best genotype differs from its called one.
     */
    bool rescore(StagedSiteTable& staged, const PhaseTable& phases, const ReadStrandTable& strands,
                 TemperFit& fit,
                 const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
                 bool show_calibration);

private:
    /// Return the calibration for this round, from the parameters, the earlier fit, or a new fit.
    ///
    /// The calibration is also stored in `counters`.
    StrandCalibration choose_calibration(const vector<PhaseSite>& phase_sites,
                                         const unordered_set<size_t>& phase_flips,
                                         const LambdaTable& lambda);

    /// Write what the last `rescore()` did, under `calibration`, to standard error.
    void report(const StrandCalibration& calibration, bool show_calibration) const;

    /// Whether re-genotyping is on.
    bool on = false;
    /// The parameters of re-genotyping.
    RegenotypeParams regenotype_params;
    /// The most linkage passes.
    size_t max_passes = 2;
    /// The ledger file, or empty for none.
    string ledger_path;
    /// Finds the genotyper that recomputes GQ.
    SiteReader reader;
    /// What the last `rescore()` did.
    RegenotypeCounters counters;
};

}

#endif
