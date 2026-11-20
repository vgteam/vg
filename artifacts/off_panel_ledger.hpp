#ifndef VG_OFF_PANEL_LEDGER_HPP_INCLUDED
#define VG_OFF_PANEL_LEDGER_HPP_INCLUDED

#include <memory>

#include "off_panel.hpp"
#include "site_read_source.hpp"

namespace vg::multipass {

using namespace std;

/**
 * The record of what became of each off-panel edit candidate, which writes the summary.
 *
 * The read source counts the candidates as it fetches reads. The direct pass files the outcome of
 * each candidate the call bar admits, from whichever thread genotyped it. Once every record is
 * written, the ledger reports the counts.
 */
class OffPanelLedger {
public:
    /// Turn off-panel detection on, for the candidates `reads` counts.
    ///
    /// `bar` decides which candidates become edit sites, and which edit sites are staged.
    void configure(const SiteReadSource* reads, const EditCallBar& bar);

    /// Return whether off-panel detection is on.
    bool is_enabled() const { return candidate_source != nullptr; }

    /// Start filing outcomes, from thread numbers below `threads`.
    void open(size_t threads);

    /// Get the read source that counts the candidates, or null when detection is off.
    const SiteReadSource* source() const { return candidate_source; }

    /// Get the call bar.
    const EditCallBar& bar() const { return call_bar; }

    /// Get the outcomes, for the direct pass to file into.
    ///
    /// @return Null before `open`.
    EditOutcomes* outcomes() { return outcomes_filed.get(); }

    /**
     * Report what became of the candidates on standard error.
     *
     * Do nothing when detection is off.
     *
     * @param nodes_after_pass The number of nodes whose top-level edit sites were genotyped after
     * the direct pass.
     */
    void report(size_t nodes_after_pass) const;

private:
    /// The read source that counts the candidates, or null when detection is off.
    const SiteReadSource* candidate_source = nullptr;
    /// Which candidates become edit sites, and which edit sites are staged.
    EditCallBar call_bar;
    /// What became of each candidate the call bar admitted.
    unique_ptr<EditOutcomes> outcomes_filed;
};

}

#endif
