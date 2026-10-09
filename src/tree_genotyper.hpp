#ifndef VG_TREE_GENOTYPER_HPP_INCLUDED
#define VG_TREE_GENOTYPER_HPP_INCLUDED

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"
#include "candidate_finder.hpp"
#include "child_placer.hpp"
#include "genotype_linker.hpp"
#include "ploidy_regions.hpp"
#include "site_genotyper.hpp"
#include "staged_site.hpp"

namespace vg {

using namespace std;

/**
 * The direct pass for one top-level site. It genotypes the site from its reads, gives it to the
 * linkage model and stages it, and then genotypes the sites below it: with nested calling, each
 * child chain its called alleles cross (see `ChildPlacer`); with --top-down, each child against
 * the traversals its parent's called alleles allow. Each site below is staged in the same way, and
 * so on down. Records are built later, from the staged sites (see `RecordRenderer`).
 *
 * Several threads may genotype different top-level sites at once.
 */
class TreeGenotyper {
public:
    /// What the genotyper reads and writes. None of it is owned.
    struct Parts {
        const PathPositionHandleGraph* graph = nullptr;
        /// The snarls, for a --top-down site's children.
        SnarlManager* snarl_manager = nullptr;
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
        /// A site's record key (see `VCFOutputCaller::record_key_of`).
        function<size_t(const Snarl&)> record_key_of;
    };

    /// Which sites below a top-level site are genotyped.
    struct Options {
        /// Nested calling: genotype the child chains the called alleles cross.
        bool nested_calling = false;
        /// With nested calling, also genotype child chains the reference does not cross, with no
        /// line.
        bool off_reference = false;
        /// --top-down: genotype each child against the traversals its parent's called alleles
        /// allow.
        bool top_down = false;
        /// -Y: under --top-down, a parent allele that skips a child gives it a star allele rather
        /// than a missing one.
        bool star_allele = false;
    };

    void configure(const Parts& parts, const Options& options);

    /// Genotype and stage the top-level site `managed_snarl` and the sites below it. Returns
    /// false when the site itself could not be genotyped, so that the walk can genotype its
    /// children as top-level sites instead.
    bool genotype(const Snarl& managed_snarl);

private:
    /// Genotype and stage one site, then the sites below it.
    /// @param parent_ref_path_name Reference path from parent (for off-reference snarls)
    /// @param parent_ref_interval Reference interval from parent
    /// @param parent_child_trav_sets If non-null, contains one TraversalSet per parent allele.
    ///                               Each set contains all traversals through this child that are
    ///                               consistent with that parent allele, and the child's genotype
    ///                               takes one allele from each set. --top-down passes them;
    ///                               nested calling passes null.
    /// @param ploidy_override If >= 0, the ploidy to genotype this snarl at, instead of the
    ///                        contig's or the --ploidy-bed region's. Nested calling passes the
    ///                        number of the parent's called alleles that cross the child, or the
    ///                        parent's ploidy when none does.
    /// @param placement Where the snarl sits in the nesting tree: the default for a top-level
    ///                  snarl, and what `ChildPlacer::place` gave a child. A --top-down child
    ///                  takes its parent's.
    bool genotype_tree(const Snarl& managed_snarl, const string& parent_ref_path_name,
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
    unique_ptr<StagedSite> stage_render_record(const Snarl& snarl,
                                               const vector<int>& trav_genotype, int ref_trav_idx,
                                               unique_ptr<SnarlCaller::CallInfo>& call_info,
                                               SiteScore* score, const string& ref_path_name,
                                               int ref_offset, int ploidy) const;

    /// Give a staged site what it keeps of the decomposition: `children`, which `snarl`, the
    /// site, holds, whether it is a leaf, and its chain.
    void fill_tree_fields(const Snarl& snarl, const SiteChildren& children,
                          StagedSite& site) const;

    Parts parts;
    Options options;
};

}

#endif
