#ifndef VG_EDIT_SITE_PLACER_HPP_INCLUDED
#define VG_EDIT_SITE_PLACER_HPP_INCLUDED

#include <functional>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include "descent_policy.hpp"
#include "edit_snarl_decomposition_overlay.hpp"
#include "../src/handle.hpp"
#include "off_panel.hpp"
#include "reference_path_table.hpp"
#include "site_read_source.hpp"
#include "site_values.hpp"

namespace vg::multipass {

using namespace std;

/// The edit site of an off-panel candidate, as a child of the site holding it.
struct HeldEdit {
    /// The edit site's net handle, as an edit in the overlay of sites.
    net_handle_t net;
    /// The alleles of the edit site.
    CandidateAlleles alleles;
    EditCandidate candidate;
};

/// Make the outcome of a candidate for the dump, saying where it was placed, or why it was not.
EditOutcome unplaced_outcome(const EditCandidate& c, const string& why);

/**
 * A placer of off-panel edit sites among the sites of the decomposition.
 *
 * It takes the candidates the call bar admits, one per base, and finds the site that holds the
 * nodes of each. It makes the candidate's edit site, and adds its edit to the overlay of sites, as
 * a child of the lowest element holding its nodes. A candidate on a node that no site holds makes a
 * free edit site, which is genotyped as a top-level site. The placer sees that each node's free
 * edit sites are placed once.
 *
 * Several threads may place edit sites at once.
 */
class EditSitePlacer {
public:
    /// The parts the placer reads. None is owned.
    struct Parts {
        const PathPositionHandleGraph* graph = nullptr;
        /// Finds the sites that hold a node.
        const DescentPolicy* descent_policy = nullptr;
        /// Off-panel detection's candidates.
        const SiteReadSource* edit_source = nullptr;
        /// What became of each candidate the call bar admits.
        EditOutcomes* edit_outcomes = nullptr;
        /// The reference paths.
        const ReferencePathTable* ref_paths = nullptr;
        /// The sites, which get each edit site's edit.
        EditSnarlDecompositionOverlay* sites = nullptr;
    };

    /// Set the parts the placer reads, and the bar a candidate must clear.
    void configure(const Parts& parts, const EditCallBar& bar) {
        this->parts = parts;
        this->bar = bar;
    }

    /**
     * List the edit sites that a site holds itself, rather than through a child site.
     *
     * Read each one's node the way the walks of the site read it, preferring the reference walk.
     *
     * @param site The site, as its net handle in the overlay of sites.
     * @param walks The candidate walks of the site.
     * @param ref_trav_idx The index in `walks` of the reference walk.
     */
    vector<HeldEdit> held_by(const net_handle_t& site, const vector<Traversal>& walks,
                             int ref_trav_idx) const;

    /**
     * Place the free edit sites of `candidates`, all on one node, unless that node's were already
     * placed.
     *
     * Add the edit of each to the overlay of sites, and call `genotype` on it.
     */
    void place_free_edits(vector<EditCandidate> candidates,
                          const function<void(const HeldEdit& edit)>& genotype);

    /**
     * Group the candidates the call bar admits, on nodes whose free edit sites are not placed yet,
     * by node.
     *
     * @param candidates Sorted by node, offset, kind and ALT.
     */
    vector<vector<EditCandidate>> group_unclaimed(vector<EditCandidate> candidates);

    /// Count the nodes `group_unclaimed` last grouped.
    size_t nodes_after_pass() const { return after_pass_nodes; }

    /// Append to `out` the nodes that no site holds along a reference path from `bound`.
    ///
    /// Go both ways, until a node a site holds, and at most a fixed number of nodes each way.
    void free_nodes_along_reference(nid_t bound, vector<nid_t>& out) const;

private:
    /// Pick the candidates the call bar admits, one per base: the ALT the most fragments carry.
    ///
    /// @param candidates Sorted by node, offset, kind and ALT.
    vector<EditCandidate> admitted(const vector<EditCandidate>& candidates) const;

    /// Return true if every node of the walk of `c` lies directly in the site `where` gives.
    bool walk_in_one_site(const EditCandidate& c, const NodePlacement& where) const;

    /// Return true if the reference reads the walk of `c` forward and once, or not at all.
    bool reference_reads_whole_walk(const EditCandidate& c) const;

    /// Add `site` to the overlay of sites as an edit, and return its net handle.
    net_handle_t add_edit(const EditSite& site) const;

    /**
     * Make the edit site of `c` into `out`, read the way `read_as` reads its node.
     *
     * @param where The place of the node of `c`.
     * @return False if the change has a walk of several nodes, or is an indel, and `read_as` is
     * backward. Also false if the walk leaves the site that holds the node, or a reference path
     * reads the walk only in part, backward, or more than once.
     */
    bool edit_site_of(const EditCandidate& c, handle_t read_as, const NodePlacement& where,
                      EditSite& out) const;

    Parts parts;
    EditCallBar bar;
    /// Guards `free_edit_nodes`.
    std::mutex free_edit_mutex;
    /// The nodes whose free edit sites have been placed.
    unordered_set<nid_t> free_edit_nodes;
    /// See `nodes_after_pass()`.
    size_t after_pass_nodes = 0;
};

}

#endif
