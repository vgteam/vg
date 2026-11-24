#ifndef VG_READ_STRAND_TABLE_HPP_INCLUDED
#define VG_READ_STRAND_TABLE_HPP_INCLUDED

#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "linkage_model.hpp"
#include "read_phasing.hpp"
#include "regenotype.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A record of the temper that the first re-genotyping round fitted, kept for the later rounds.
 *
 * The temper describes how reliable the summed strand log-odds of the reads are, which does not
 * depend on the genotypes called, so it is fitted once. The record also holds the ceiling the
 * temper was fitted against, and the calibration table of the fit.
 */
struct TemperFit {
    /// The calibration the first round used.
    StrandCalibration calibration;
    /// The |Lambda| of each bin of the calibration table.
    vector<double> abs_lambda;
    /// The observed agreement with the chain in each bin.
    vector<double> observed;
    /// The predicted agreement with the chain in each bin.
    vector<double> predicted;
    /// The number of reads in each bin.
    vector<size_t> count;

    /// Keep the fit that a round's counters hold.
    void keep(const RegenotypeCounters& counters);

    /// Put the kept fit into a round's counters.
    void restore(RegenotypeCounters& counters) const;
};

/**
 * A table of which strand each read came from, as read phasing last decided the sites.
 *
 * Read phasing fills `sites()` with the evidence of each read at each diploid heterozygous site,
 * reduced to the site's chosen pair of alleles, and fills `flips()` with the sites whose pair it
 * reversed. Re-genotyping reads both. For the anchors, `build_lambda` then sums the strand
 * log-odds of each read over the sites, and `read_strand_log_odds` answers for one read at one
 * site.
 */
class ReadStrandTable {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// Get the evidence of each diploid heterozygous site, reduced to its chosen pair.
    vector<PhaseSite>& sites() { return phase_sites; }
    /// Get the evidence of each diploid heterozygous site, reduced to its chosen pair.
    const vector<PhaseSite>& sites() const { return phase_sites; }

    /// Get the record keys of the sites whose chosen pair read phasing reversed.
    ///
    /// The contribution of such a site to the strand log-odds of a read has the opposite sign.
    unordered_set<size_t>& flips() { return phase_flips; }
    /// Get the record keys of the sites whose chosen pair read phasing reversed.
    const unordered_set<size_t>& flips() const { return phase_flips; }

    /// Get the sites that the reads ordered where the panel did not (see `PhaseSite::panel_ordered`).
    unordered_set<size_t>& read_ordered() { return ordered_by_reads; }
    /// Get the sites that the reads ordered where the panel did not (see `PhaseSite::panel_ordered`).
    const unordered_set<size_t>& read_ordered() const { return ordered_by_reads; }

    /**
     * Sum the strand log-odds of each read over `sites()`, for the anchors.
     *
     * Call this once `sites()` and `flips()` are final. The table is built here, and not taken
     * from re-genotyping, because re-genotyping may not have run, and builds its table before the
     * sites are final.
     *
     * @param calls The phase set of each site, as `PhaseCall`s, of which the last one for a
     * record wins.
     * @param phase_set_id Numbers a phase set, by its contig and its phase set in that contig.
     * @param fitted The calibration re-genotyping used, or one with a temper of 0 or less to fit
     * one here under `params`.
     */
    void build_lambda(const vector<PhaseCall>& calls,
                      const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
                      const StrandCalibration& fitted, const RegenotypeParams& params);

    /**
     * Return the tempered strand log-odds of a read, leaving out the site `record_key`.
     *
     * Leaving out the site keeps a site from judging its own reads. A positive value favors anchor
     * slot 0. The value is 0 when there is no table, no fitted temper, or no other site the read
     * contributed to. It is NaN when the read has a strand that is not usable in the phase set of
     * the site (see `read_strand_usable`); a split homozygous site drops such a read, but places a
     * read with no strand by a coin.
     *
     * @param site_own The log-odds of each read at this site alone, from
     * `site_own_strand_log_odds` for the same `record_key`, or null to compute them here. A caller
     * looking up many reads at one site can build them once.
     */
    double read_strand_log_odds(size_t record_key, std::string_view read_name,
                                const unordered_map<uint64_t, double>* site_own = nullptr) const;

    /**
     * Fill `out` with the log-odds of each read at the site `record_key` alone.
     *
     * These are what `read_strand_log_odds` leaves out.
     *
     * @return False, leaving `out` alone, when there is no table, no fitted temper, or the site
     * contributed nothing.
     */
    bool site_own_strand_log_odds(size_t record_key, unordered_map<uint64_t, double>& out) const;

private:
    /// The evidence of each diploid heterozygous site, reduced to its chosen pair.
    vector<PhaseSite> phase_sites;
    /// The record keys of the sites whose chosen pair read phasing reversed.
    unordered_set<size_t> phase_flips;
    /// The record keys of the sites that the reads ordered where the panel did not.
    unordered_set<size_t> ordered_by_reads;

    /// The strand log-odds of each read, from `build_lambda`, by hash of the read name.
    ///
    /// Positive means strand 0 of the phase set of the read, which is GT field 0 and anchor slot
    /// 0. The key is the read alone, because a homozygous site, which has no PhaseSite, looks up
    /// reads in it too.
    LambdaTable lambda;
    /// The entry of `phase_sites` for each record key.
    ///
    /// These point into `phase_sites`, which must not change after `build_lambda`.
    unordered_map<size_t, const PhaseSite*> lambda_site;
    /// The phase set of each site, by record key.
    ///
    /// The strand of a read is usable only at sites of the phase set it was found in (see
    /// `read_strand_usable`).
    unordered_map<size_t, size_t> lambda_phase_set;
    /// The calibration of the strand log-odds, with temper 0 when no fit was possible.
    ///
    /// At temper 0, the tempered strand log-odds of every read are 0, so no read counts as placed.
    StrandCalibration calibration;
};

}

#endif
