#ifndef VG_OFF_PANEL_LEDGER_HPP_INCLUDED
#define VG_OFF_PANEL_LEDGER_HPP_INCLUDED

#include <memory>
#include <string>

#include "../src/handle.hpp"
#include "off_panel.hpp"
#include "reference_path_table.hpp"
#include "site_read_source.hpp"
#include "tree_genotyper.hpp"

namespace vg::multipass {

using namespace std;

/**
 * The record of what became of each off-panel edit candidate, which writes the dump and the
 * summary.
 *
 * The read source counts the candidates as it fetches reads. The direct pass files the outcome of
 * each candidate the call bar admits, from whichever thread genotyped it. Once every record is
 * written, the ledger reports the counts, and writes each candidate with its outcome to the dump,
 * placed on the reference.
 */
class OffPanelLedger {
public:
    /// Make a ledger that places candidates on `graph`, on the reference paths of `ref_paths`,
    /// shifted by their offsets.
    ///
    /// Neither is owned.
    OffPanelLedger(const PathPositionHandleGraph& graph, const ReferencePathTable& ref_paths);

    /// Turn off-panel detection on, for the candidates `reads` counts.
    ///
    /// `bar` decides which candidates become edit sites, and which edit sites are staged. With a
    /// `dump` path, every candidate is written there at the end.
    void configure(const SiteReadSource* reads, string dump, const EditCallBar& bar);

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
     * Report what became of the candidates, and write every candidate to the dump, if there is
     * one.
     *
     * Do nothing when detection is off.
     *
     * @param nodes_after_pass The number of nodes whose top-level edit sites were genotyped after
     * the direct pass, for the report.
     */
    void write(size_t nodes_after_pass) const;

private:
    /// Report what became of the candidates.
    void report(size_t nodes_after_pass) const;

    /// Find where a node's base lies on the reference.
    ///
    /// @return The contig and the 1-based position of the base at forward offset `offset`, or an
    /// empty contig when no reference path visits the node.
    pair<string, int64_t> reference_position(nid_t node, uint32_t offset) const;

    const PathPositionHandleGraph& graph;
    /// The reference paths the candidates are placed on.
    const ReferencePathTable& ref_paths;

    /// The read source that counts the candidates, or null when detection is off.
    const SiteReadSource* candidate_source = nullptr;
    /// Where to write the dump, or empty for nowhere.
    string dump_path;
    /// Which candidates become edit sites, and which edit sites are staged.
    EditCallBar call_bar;
    /// What became of each candidate the call bar admitted.
    unique_ptr<EditOutcomes> outcomes_filed;
};

}

#endif
