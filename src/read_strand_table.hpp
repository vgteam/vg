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

namespace vg {

using namespace std;

/// The temper the first re-genotyping round fitted, and that fit's calibration table, kept for
/// the later rounds. The temper describes how reliable the reads' summed strand log-odds are, not
/// which genotypes are called, so it is fitted once.
struct TemperFit {
    double temper = 0.0;
    vector<double> abs_lambda, observed, predicted;
    vector<size_t> count;

    /// Keep the fit a round's counters hold.
    void keep(const RegenotypeCounters& counters);

    /// Start a round's counters from the kept fit. The ceiling is not part of the fit, so the
    /// counters keep their default ceiling.
    void restore(RegenotypeCounters& counters) const;
};

/**
 * What the reads say about each site's strands, as read phasing last decided them.
 *
 * Read phasing reduces each diploid heterozygous site's per-read evidence to its chosen pair
 * (`sites()`), and decides at which sites to reverse the pair (`flips()`). Re-genotyping reads
 * both. For the anchors, `build_lambda` then sums each read's strand log-odds over the sites,
 * and `read_strand_log_odds` answers for one read at one site.
 */
class ReadStrandTable {
public:
    using PhaseCall = LinkageCollector::PhaseCall;

    /// Each diploid heterozygous site's per-read evidence, reduced to its chosen pair.
    vector<PhaseSite>& sites() { return phase_sites; }
    const vector<PhaseSite>& sites() const { return phase_sites; }

    /// The record keys whose chosen pair read phasing reversed. Their sites' contributions to a
    /// read's strand log-odds enter with the opposite sign.
    unordered_set<size_t>& flips() { return phase_flips; }
    const unordered_set<size_t>& flips() const { return phase_flips; }

    /// Sum each read's strand log-odds over `sites()`, for the anchors, once `sites()` and
    /// `flips()` are final. `calls` gives each site's phase set, which `phase_set_id` numbers.
    /// The log-odds are tempered with `fitted_temper` and `fitted_ceiling` where a temper was
    /// fitted (`fitted_temper` above 0), and otherwise with a temper fitted here, under
    /// `params`.
    void build_lambda(const vector<PhaseCall>& calls,
                      const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
                      double fitted_temper, double fitted_ceiling, const RegenotypeParams& params);

    /// This read's tempered strand log-odds, leaving out `record_key`, so that a site does not judge
    /// its own reads. Positive names slot 0. Zero means none: no table, no other contributing site,
    /// or no fitted temper. NaN means the read has a strand that is not usable in the site's phase
    /// set (see `read_strand_usable`). Re-genotyping does not use this, and gives such a read 0.
    ///
    /// `site_own`, if given, is the site's own log-odds per read, from `site_own_strand_log_odds`
    /// for the same `record_key`; a caller looking up many reads at one site builds it once.
    double read_strand_log_odds(size_t record_key, std::string_view read_name,
                                const unordered_map<uint64_t, double>* site_own = nullptr) const;

    /// Fill `out` with each read's log-odds from the site `record_key` alone, which
    /// `read_strand_log_odds` leaves out. Returns false, leaving `out` alone, when there is nothing to
    /// leave out: no table, no fitted temper, or a site that contributed nothing.
    bool site_own_strand_log_odds(size_t record_key, unordered_map<uint64_t, double>& out) const;

protected:
    vector<PhaseSite> phase_sites;
    unordered_set<size_t> phase_flips;

    /// Each read's strand log-odds, from `build_lambda`. Positive means strand 0 of the read's
    /// phase set, which is GT field 0 and anchor slot 0. Keyed by read alone, as a homozygous
    /// site, which has no PhaseSite, needs.
    LambdaTable lambda;
    /// `phase_sites` indexed by record, so that a site's own contribution can be subtracted.
    /// Points into `phase_sites`, which must not be rebuilt afterwards.
    unordered_map<size_t, const PhaseSite*> lambda_site;
    /// Each site's phase set, by record, since a read's strand is usable only at sites of the
    /// phase set it was found in (see `read_strand_usable`).
    unordered_map<size_t, size_t> lambda_phase_set;
    /// The temper. Zero, when no fit was possible, makes every read's tempered strand log-odds
    /// zero, so no read counts as placed.
    double lambda_temper = 0.0;
    double lambda_ceiling = 1.0;
};

}

#endif
