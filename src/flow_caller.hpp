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
#include "candidate_finder.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "graph_caller.hpp"
#include "vcf_output_caller.hpp"
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
class FlowCaller : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
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

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

    /// Do not genotype a snarl with more edges than this, including those of nested snarls
    /// (--max-snarl-edges). `call_top_level_snarls` then genotypes the snarl's children as if they
    /// were top-level snarls. Zero removes the limit, which is also the default.
    void set_max_snarl_edges(size_t edges) { candidates.set_max_snarl_edges(edges); }

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

    /// Finds each site's reference path and candidate traversals. Declared after the members it
    /// reads.
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
                             const ChildTraversalSets* parent_child_trav_sets = nullptr);

};

}

#endif
