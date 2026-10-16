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

#include "handle.hpp"
#include "snarls.hpp"
#include "candidate_finder.hpp"
#include "child_placer.hpp"
#include "genotype_linker.hpp"
#include "off_panel.hpp"
#include "site_read_source.hpp"
#include "ploidy_regions.hpp"
#include "site_genotyper.hpp"
#include "staged_site.hpp"

namespace vg {
namespace multipass {

using namespace std;

/**
 * The direct pass for one top-level site. It genotypes the site from its reads, gives it to the
 * linkage model and stages it, and then genotypes the sites below it: with nested calling, each
 * child chain its called alleles cross (see `ChildPlacer`); with top-down calling, each child
 * against the traversals its parent's called alleles allow. Each site below is staged in the same
 * way, and so on down. Records are built later, from the staged sites (see `RecordRenderer`).
 *
 * Several threads may genotype different top-level sites at once.
 */
class TreeGenotyper {
public:
    /// What the genotyper reads and writes. None of it is owned.
    struct Parts {
        const PathPositionHandleGraph* graph = nullptr;
        const CandidateFinder* candidates = nullptr;
        const SiteGenotyper* genotyper = nullptr;
        const GenotypeLinker* linker = nullptr;
        StagedSiteTable* staged_sites = nullptr;
        const ChildPlacer* child_placer = nullptr;
        DescentCounters* descent_counters = nullptr;
        const PloidyRegions* ploidy_regions = nullptr;
        /// The offset and ploidy of each reference path.
        const map<string, size_t>* ref_offsets = nullptr;
        const map<string, int>* ref_ploidies = nullptr;
        /// The ID of the site with bounds `site`, as its records name it (see `StagedSite::id`).
        function<string(const SiteBounds& site)> site_id;
        /// Off-panel detection's candidates, or null when it is off.
        const SiteReadSource* edit_source = nullptr;
        /// The source the sites' reads come from. A top-level site holds its reads there for the
        /// sites below it (see `SiteReadSource::hold_reads`). Null for none.
        const SiteReadSource* read_source = nullptr;
        /// What became of each candidate the call bar admits. Set when `edit_source` is.
        EditOutcomes* edit_outcomes = nullptr;
    };

    /// Which sites below a top-level site are genotyped.
    struct Options {
        /// Nested calling: genotype the child chains the called alleles cross.
        bool nested_calling = false;
        /// With nested calling, also genotype child chains the reference does not cross, with no
        /// line.
        bool off_reference = false;
        /// Top-down calling: genotype each child against the traversals its parent's called
        /// alleles allow.
        bool top_down = false;
        /// With `top_down`, a parent allele that skips a child gives it a star allele rather than
        /// a missing one.
        bool star_allele = false;
        /// Which off-panel candidates become edit sites, and which edit sites are staged.
        EditCallBar edit_bar;
    };

    void configure(const Parts& parts, const Options& options);

    /// Genotype and stage the top-level site `site` and the sites below it. Returns false when
    /// the site itself could not be genotyped, so that the walk can genotype its children as
    /// top-level sites instead.
    bool genotype(const SiteView& site);

    /// Genotype and stage, as top-level sites, the edit sites of the off-panel candidates
    /// `candidates` (sorted) that no site holds and that the direct pass did not reach: those
    /// on a node that bounds no top-level site. A top-level site genotypes those on its two
    /// boundary nodes, while their reads are in memory, and an edit site a site holds is one of
    /// that site's children, genotyped and staged in its descent. Runs on several threads.
    void genotype_top_level_edits(vector<EditCandidate> candidates);

    /// How many nodes' top-level edit sites `genotype_top_level_edits` was left to genotype.
    size_t edit_nodes_after_pass() const { return after_pass_nodes; }

private:
    /// An off-panel candidate's edit site, as a child of the site holding it.
    struct HeldEdit {
        ChildSite child;
        EditSite site;
        EditCandidate candidate;
    };

    /// Genotype and stage, as top-level sites, the edit sites of `candidates` (sorted) that no
    /// site holds, unless their node's were already, and note their node as done.
    void genotype_free_edits(vector<EditCandidate> candidates);

    /// The nodes whose free edit sites have been genotyped.
    std::mutex free_edit_mutex;
    unordered_set<nid_t> free_edit_nodes;
    size_t after_pass_nodes = 0;

    /// Append to `out` the nodes no site holds that run along a reference path from `bound`, in
    /// either direction, until a node a site holds, at most a fixed number each way.
    void free_nodes_along_reference(nid_t bound, vector<nid_t>& out) const;

    /// The candidates the call bar admits, one per base: the ALT the most fragments carry.
    vector<EditCandidate> admitted(vector<EditCandidate> candidates) const;

    /// The edit site of `c`, read the way `read_as` reads its node, into `out`. False for an
    /// indel on a node read backward, for one whose walk leaves the site `where` holds its anchor
    /// in, and for one whose walk a reference path reads only in part, backward or more than once.
    bool edit_site_of(const EditCandidate& c, handle_t read_as, const NodePlacement& where,
                      EditSite& out) const;

    /// The edit sites that the site `view`, genotyped as `bounds` with candidate walks `walks`,
    /// holds itself rather than through a child site, each read as the walks read its node,
    /// the reference walk `ref_trav_idx` first.
    vector<HeldEdit> edits_held_by(const SiteView& view, const SiteBounds& bounds,
                                   const vector<Traversal>& walks, int ref_trav_idx) const;

    /// Genotype one edit site and stage it if it is called with its ALT: as a top-level site
    /// where `top_level`, and otherwise as a child placed at `placement` under a site on the
    /// reference path `parent_ref_path_name` over `parent_ref_interval`, at `ploidy_override`.
    void genotype_edit_site(const HeldEdit& edit, const vector<SiteBounds>& enclosing,
                            const string& parent_ref_path_name,
                            pair<size_t, size_t> parent_ref_interval, int ploidy_override,
                            const NestingPlacement& placement, bool top_level);

    /// Genotype and stage one site, then the sites below it.
    /// @param parent_ref_path_name Reference path from parent (for off-reference snarls)
    /// @param parent_ref_interval Reference interval from parent
    /// @param parent_child_trav_sets If non-null, contains one TraversalSet per parent allele.
    ///                               Each set contains all traversals through this child that are
    ///                               consistent with that parent allele, and the child's genotype
    ///                               takes one allele from each set. Top-down calling passes
    ///                               them; nested calling passes null.
    /// @param ploidy_override If >= 0, the ploidy to genotype this snarl at, instead of the
    ///                        contig's or a per-region override's. Nested calling passes the
    ///                        number of the parent's called alleles that cross the child, or the
    ///                        parent's ploidy when none does.
    /// @param placement Where the snarl sits in the nesting tree: the default for a top-level
    ///                  snarl, and what `ChildPlacer::place` gave a child. A child in
    ///                  top-down calling takes its parent's.
    bool genotype_tree(const SiteView& view, const string& parent_ref_path_name,
                       pair<size_t, size_t> parent_ref_interval,
                       const ChildTraversalSets* parent_child_trav_sets, int ploidy_override,
                       const NestingPlacement& placement);

    /// Genotype one site at `ploidies`. `score` is set to the score inside the returned call
    /// info.
    pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> genotype_site(
        const SiteBounds& site, const vector<Traversal>& travs, int ref_trav_idx,
        const Ploidies& ploidies, const vector<SiteBounds>& enclosing,
        const string& ref_path_name, pair<size_t, size_t> ref_range, SiteScore*& score) const;

    /// Make a site's `StagedSite` from its genotype, moving `call_info` into it. `travs` is left
    /// empty, because the sites below still read the traversals; the caller moves them in once
    /// those are done.
    unique_ptr<StagedSite> stage_render_record(const SiteBounds& bounds, const string& site_id,
                                               const vector<int>& trav_genotype, int ref_trav_idx,
                                               unique_ptr<SnarlCaller::CallInfo>& call_info,
                                               SiteScore* score, const string& ref_path_name,
                                               int ref_offset, int ploidy) const;

    /// Give a staged site what it keeps of the decomposition: `children`, which `view`, the
    /// site, holds, whether it is a leaf, its chain, and the sites enclosing it.
    void fill_tree_fields(const SiteView& view, const SiteChildren& children,
                          StagedSite& site) const;

    Parts parts;
    Options options;
};

}
}

#endif
