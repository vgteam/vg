#include <algorithm>
#include <iostream>

#include "off_panel_ledger.hpp"

namespace vg::multipass {

void OffPanelLedger::configure(const SiteReadSource* reads, const EditCallBar& bar) {
    candidate_source = reads;
    call_bar = bar;
}

void OffPanelLedger::open(size_t threads) {
    outcomes_filed = std::make_unique<EditOutcomes>(threads);
}

void OffPanelLedger::report(size_t nodes_after_pass) const {
    if (!is_enabled()) {
        return;
    }
    size_t top = 0;
    size_t nested = 0;
    size_t staged = 0;
    size_t unplaced = 0;
    size_t indels_staged = 0;
    for (const EditOutcome& o : outcomes_filed->all()) {
        if (o.placed == "top") {
            ++top;
        } else if (o.placed == "nested") {
            ++nested;
        } else {
            ++unplaced;
        }
        if (o.staged) {
            ++staged;
            if (o.kind == 'D' || o.kind == 'I') {
                ++indels_staged;
            }
        }
    }
    const vector<EditCandidate>& candidates = candidate_source->counted_edit_candidates();
    cerr << "[vg call] off-panel: " << candidates.size()
         << " candidates counted; the call bar's genotyped as " << top << " top-level and "
         << nested << " nested sites, " << unplaced << " left unplaced; " << staged
         << " staged with their ALT; " << nodes_after_pass
         << " nodes' top-level sites left to after the direct pass" << endl;
    const auto indels = std::ranges::count_if(candidates,
                                              [](const EditCandidate& c) { return c.is_indel(); });
    cerr << "[vg call] off-panel: " << indels << " of the candidates are indels; "
         << indels_staged << " indels staged with their ALT" << endl;
}

}
