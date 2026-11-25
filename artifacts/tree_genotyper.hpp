#ifndef VG_TREE_GENOTYPER_HPP_INCLUDED
#define VG_TREE_GENOTYPER_HPP_INCLUDED

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <string>
#include <utility>
#include <vector>

#include "../src/handle.hpp"
#include "../src/snarls.hpp"
#include "candidate_finder.hpp"
#include "descent_policy.hpp"
#include "edit_site_placer.hpp"
#include "edit_snarl_decomposition_overlay.hpp"
#include "genotype_linker.hpp"
#include "off_panel.hpp"
#include "site_read_source.hpp"
#include "ploidy_regions.hpp"
#include "reference_path_table.hpp"
#include "read_likelihood_caller.hpp"
#include "staged_site.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A genotyper that runs the direct pass on one top-level site and the sites below it.
 *
 * It genotypes the site from its reads, files it in the linkage model and stages it. Then it
 * genotypes the sites below it in the same way, and so on down. With nested calling, those are
 * the child chains that the called alleles cross (see `DescentPolicy`).
 *
 * Several threads may genotype different top-level sites at once.
 */
class TreeGenotyper {
public:
    /// The objects the genotyper reads and writes. None of them is owned.
    struct Parts {
        const PathPositionHandleGraph* graph = nullptr;
        /// Finds the reference path and candidate traversals of each site.
        const CandidateFinder* candidates = nullptr;
        /// Genotypes each site from its reads.
        ReadLikelihoodSnarlCaller* genotyper = nullptr;
        /// Files each site in the linkage model.
        const GenotypeLinker* linker = nullptr;
        /// Holds the staged sites.
        StagedSiteTable* staged_sites = nullptr;
        /// Lists the child chains to genotype under a site.
        const DescentPolicy* descent_policy = nullptr;
        DescentCounters* descent_counters = nullptr;
        const PloidyRegions* ploidy_regions = nullptr;
        /// The reference paths, with their offsets and ploidies.
        const ReferencePathTable* ref_paths = nullptr;
        /// Makes the ID of the site with start bound `start` and end bound `end`, as its records
        /// name it (see `StagedSite::id`).
        function<string(const handle_t& start, const handle_t& end)> site_id;
        /// Off-panel detection's candidates, or null when it is off.
        const SiteReadSource* edit_source = nullptr;
        /// The source the reads of the sites come from, or null for none.
        ///
        /// A top-level site holds its reads there for the sites below it (see
        /// `SiteReadSource::hold_reads`).
        const SiteReadSource* read_source = nullptr;
        /// What became of each candidate the call bar admits. Set when `edit_source` is.
        EditOutcomes* edit_outcomes = nullptr;
        /// The sites, with the edit sites laid over them. Each edit site is added there when it
        /// is placed, and stays for the rest of the run.
        EditSnarlDecompositionOverlay* sites = nullptr;
    };

    /// The choices of which sites below a top-level site are genotyped.
    struct Options {
        /// Nested calling: genotype the child chains the called alleles cross.
        bool nested_calling = false;
        /// With nested calling, also genotype child chains the reference does not cross, with no
        /// line.
        bool off_reference = false;
        /// Whether groups of each site's walks, from `candidate_groups`, choose its genotype. The
        /// walks then choose within the chosen groups.
        bool group_candidates = false;
        /// Which off-panel candidates become edit sites, and which edit sites are staged.
        EditCallBar edit_bar;
    };

    /// Give the genotyper the objects it works with and its choices.
    void configure(const Parts& parts, const Options& options);

    /**
     * Genotype and stage the top-level site `site` and the sites below it.
     *
     * With off-panel detection, also genotype the edit sites that no site holds on the boundary
     * nodes of the site and on the nodes along the reference beside it, while their reads are in
     * memory.
     *
     * @return False if the site itself could not be genotyped, so that the walk can genotype its
     * children as top-level sites instead.
     */
    bool genotype(const net_handle_t& site);

    /**
     * Genotype and stage the free edit sites that the direct pass did not reach, as top-level
     * sites.
     *
     * A free edit site is the edit site of an off-panel candidate on a node that no site holds.
     * The direct pass genotypes those near a top-level site (see `genotype`). This gets the rest.
     * It runs on several threads.
     *
     * @param candidates The off-panel candidates, sorted.
     */
    void genotype_top_level_edits(vector<EditCandidate> candidates);

    /// Count the nodes whose free edit sites `genotype_top_level_edits` genotyped.
    size_t edit_nodes_after_pass() const { return placer.nodes_after_pass(); }

private:
    struct SiteCall;

    /// Genotype and stage, as top-level sites, the free edit sites of `candidates`, all on one
    /// node, unless that node's were already done.
    ///
    /// @param candidates The off-panel candidates on one node, sorted.
    void genotype_free_edits(vector<EditCandidate> candidates);

    /// Find the reference path through `node`, preferring `preferred`, and the position of the
    /// node on it, into `node_start`.
    ///
    /// @return The path's name, or empty if no reference path visits the node.
    string reference_through(handle_t node, const string& preferred, size_t& node_start) const;

    /**
     * Genotype one edit site, and stage it if it is called with its ALT.
     *
     * @param parent_ref_path_name The reference path of the parent, for a nested edit site.
     * @param parent_ref_interval The interval of the parent on its reference path.
     * @param ploidy_override If not negative, the ploidy to genotype at.
     * @param placement The place of the edit site in the nesting tree.
     * @param top_level Stage the site as a top-level site, rather than as a child.
     */
    void genotype_edit_site(const HeldEdit& edit, const string& parent_ref_path_name,
                            pair<size_t, size_t> parent_ref_interval, int ploidy_override,
                            const NestingPlacement& placement, bool top_level);

    /**
     * Genotype and stage one site, then the sites below it.
     *
     * @param parent_ref_path_name The reference path of the parent, which a site off the
     * reference takes.
     * @param parent_ref_interval The interval of the parent on its reference path.
     * @param ploidy_override If not negative, the ploidy to genotype at, in place of the
     * contig's or the region's. Nested calling passes the number of the parent's called alleles
     * that cross the child, or the parent's ploidy when none does.
     * @param placement The place of the site in the nesting tree. It is the default for a
     * top-level site.
     * @return False if the site could not be genotyped.
     */
    bool genotype_tree(const net_handle_t& view, const string& parent_ref_path_name,
                       pair<size_t, size_t> parent_ref_interval, int ploidy_override,
                       const NestingPlacement& placement);

    /// Genotype a nested chain at the ploidy its parent gives it, and record and stage it as its
    /// placement says.
    ///
    /// @return False if the site could not be genotyped.
    bool genotype_nested_chain(SiteCall& call);

    /// File the staged site of a nested chain in the linkage model, unless no called parent
    /// allele crosses it.
    void record_nested_chain(SiteCall& call) const;

    /// Genotype a top-level site at its region's ploidy, and record and stage it.
    ///
    /// @return False if the site could not be genotyped.
    bool genotype_top_level_site(SiteCall& call);

    /// Stage a site to go to a render queue, and file it in the linkage model.
    void record_and_stage(SiteCall& call, unique_ptr<SnarlCaller::CallInfo>& trav_call_info,
                          SiteScore* trav_score) const;

    /// Genotype the child chains that the called alleles of a site cross, and its held edit
    /// sites.
    void descend_nested(const SiteCall& call);

    /// Merge the edit sites a site holds itself into the list of its child sites.
    ///
    /// Put each edit site before the first child site that starts after the node of the edit, so
    /// that the children are visited along the node IDs. Then the reads of each window are still
    /// in memory when its edit sites are genotyped.
    ///
    /// @pre `held` is sorted by node.
    DescentPolicy::Nested with_held_edits(const DescentPolicy::Nested& nested,
                                        const vector<HeldEdit>& held) const;

    /// Genotype one child that `DescentPolicy::place` placed under the site of `call`.
    ///
    /// A child inside a node is one of the edit sites in `held`. Mark it in `edit_placed`.
    void genotype_placed_child(const SiteCall& call, const DescentPolicy::Placed& child,
                               const vector<HeldEdit>& held, vector<bool>& edit_placed);

    /// Give the staged site of `call`, if any, its tree fields, and add it to the staged sites.
    ///
    /// @pre The sites below have been genotyped, since they read the traversals of the staged
    /// site.
    void file_staged(SiteCall& call);

    /// Genotype the site of `call` at `ploidies`, choosing among `SiteCall::walks`.
    ///
    /// @param score Set to the score inside the returned call info.
    /// @param group_of If not null, the group of each of the walks. The groups then choose the
    /// genotype, and the walks choose within the chosen groups.
    pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> genotype_site(
        const SiteCall& call, const Ploidies& ploidies, SiteScore*& score,
        const vector<int>* group_of = nullptr) const;

    /// Return the `candidate_groups` of the walks of `call` if `Options::group_candidates` is set,
    /// and nothing otherwise.
    vector<int> group_walks(const SiteCall& call) const;

    /// Make the `StagedSite` of a site from its walks and genotype, moving `call_info` into it.
    ///
    /// The staged site is at the top level of the nesting tree until it is given a placement.
    unique_ptr<StagedSite> stage_render_record(const net_handle_t& net, const string& site_id,
                                               const vector<int>& trav_genotype,
                                               vector<Traversal> travs, int ref_trav_idx,
                                               unique_ptr<SnarlCaller::CallInfo>& call_info,
                                               SiteScore* score,
                                               const ReferencePlacement& reference,
                                               int ploidy) const;

    Parts parts;
    Options options;
    /// Places the off-panel edit sites.
    EditSitePlacer placer;
};

}

#endif
