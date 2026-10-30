#ifndef VG_FLOW_CALLER_HPP_INCLUDED
#define VG_FLOW_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "../artifacts/candidate_finder.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "graph_caller.hpp"
#include "site_emitter.hpp"
#include "sorted_vcf_writer.hpp"
#include "gaf_output_caller.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/**
 * FlowCaller : Uses any traversals finder (ex, FlowTraversalFinder) to find
 * traversals, and calls those based on how much support they have.
 * Should work on any graph but will not
 * report cyclic traversals.  Supports nested calling when enabled with the
 * nested flag, recursively processing child snarls.
 * Designed to replace LegacyCaller, as it should miss fewer obviously
 * good traversals, and is not dependent on old protobuf-based structures.
 */
class FlowCaller : public GraphCaller, public GAFOutputCaller {
public:
    /// Original constructor for non-nested mode
    FlowCaller(const PathPositionHandleGraph& graph,
               SupportBasedSnarlCaller& snarl_caller,
               SnarlManager& snarl_manager,
               const string& sample_name,
               TraversalFinder& traversal_finder,
               const vector<string>& ref_paths,
               const vector<size_t>& ref_path_offsets,
               const vector<int>& ref_path_ploidies,
               AlignmentEmitter* aln_emitter,
               bool traversals_only,
               bool gaf_output,
               size_t trav_padding,
               bool genotype_snarls,
               const pair<size_t, size_t>& allele_length_range);

    /// Extended constructor for nested mode with star alleles
    FlowCaller(const PathPositionHandleGraph& graph,
               SupportBasedSnarlCaller& snarl_caller,
               SnarlManager& snarl_manager,
               const string& sample_name,
               TraversalFinder& traversal_finder,
               const vector<string>& ref_paths,
               const vector<size_t>& ref_path_offsets,
               const vector<int>& ref_path_ploidies,
               AlignmentEmitter* aln_emitter,
               bool traversals_only,
               bool gaf_output,
               size_t trav_padding,
               bool genotype_snarls,
               const pair<size_t, size_t>& allele_length_range,
               bool nested,
               bool star_allele);

    virtual ~FlowCaller();

    bool call_snarl(const Snarl& snarl) override;

    /// Skip snarls with more edges than this, counting those of nested snarls.
    ///
    /// `call_top_level_snarls` then calls the children of a skipped snarl as if they were
    /// top-level snarls. Zero, the default, removes the limit.
    void set_max_snarl_edges(size_t edges);

    /// The sorted VCF the records go to.
    SortedVCFWriter vcf_writer;

    /// Writes the record of each called site into `vcf_writer`.
    SiteEmitter site_emitter;

protected:

    /// the graph
    const PathPositionHandleGraph& graph;

    /// the traversal finder
    TraversalFinder& traversal_finder;

    /// keep track of the reference paths
    vector<string> ref_paths;
    unordered_set<string> ref_path_set;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;
    
    /// keep traco of the ploidies (todo: just one map for all path stuff!!)
    map<string, int> ref_ploidies;

    /// alignment emitter. if not null, traversals will be output here and
    /// no genotyping will be done
    AlignmentEmitter* alignment_emitter;

    /// toggle whether to genotype or just output the traversals
    bool traversals_only;

    /// toggle whether to output vcf or gaf
    bool gaf_output;

    /// toggle whether to genotype every snarl
    /// (by default, uncalled snarls are skipped, and coordinates are flattened
    ///  out to minimize variant size -- this turns all that off)
    bool genotype_snarls;

    /// clamp calling to alleles of a given length range
    /// more specifically, a snarl is only called if
    /// 1) its largest allele is >= allele_length_range.first and
    /// 2) all alleles are < allele_length_range.second
    pair<size_t, size_t> allele_length_range;

    /// --- Nested mode members ---

    /// enable recursive calling of child snarls
    bool nested = false;

    /// use * alleles for spanning haplotypes that don't traverse nested sites
    bool star_allele = false;

    /// Finds the reference path and candidate traversals of each site. It is declared after the
    /// members it reads, so that they are built first.
    CandidateFinder candidates;

    /// Internal implementation of call_snarl that accepts parent context for nested mode
    /// When nested=true, this recursively calls children after processing the current snarl
    /// @param parent_ref_path_name Reference path from parent (for off-reference snarls)
    /// @param parent_ref_interval Reference interval from parent
    /// @param parent_child_trav_sets If non-null, contains one TraversalSet per parent allele.
    ///                               Each set contains all traversals through this child that are
    ///                               consistent with that parent allele. The child genotypes by
    ///                               picking the best pair (one from each set) based on read support.
    bool call_snarl_internal(const Snarl& snarl,
                             const string& parent_ref_path_name,
                             pair<size_t, size_t> parent_ref_interval,
                             const vector<TraversalSet>* parent_child_trav_sets = nullptr);

    /**
     * Genotype a nested site from the traversal sets that its parent's genotype gives it.
     *
     * We add the traversals of the sets to the candidates of the site, and genotype the strands
     * whose set is not empty. The other strands get a star allele, or a missing one without
     * `star_allele`.
     *
     * @param parent_sets One traversal set per strand of the parent's genotype.
     * @param call_info Set to the snarl caller's information about the call, if it was asked for
     * one.
     * @return One allele per strand.
     */
    vector<int> genotype_from_parent_sets(CandidateFinder::Site& site,
                                          const vector<TraversalSet>& parent_sets,
                                          unique_ptr<SnarlCaller::CallInfo>& call_info);

    /// Write the record of a genotyped site, as VCF or as GAF.
    ///
    /// @return false only if the VCF writer refused the record.
    bool emit_site(const CandidateFinder::Site& site, const vector<int>& genotype,
                   const unique_ptr<SnarlCaller::CallInfo>& call_info, int ploidy);

    /// Call each child of a genotyped site that an allele of the genotype visits, from the
    /// traversals of the alleles.
    void call_children(const CandidateFinder::Site& site, const vector<int>& genotype);

    /// Return, for each strand of the genotype of a site, the traversals of `child` that agree
    /// with the traversal of its allele. A star or missing allele gets an empty set.
    vector<TraversalSet> child_traversal_sets(const CandidateFinder::Site& site,
                                              const vector<int>& genotype,
                                              const Snarl& child) const;
};

}

#endif
