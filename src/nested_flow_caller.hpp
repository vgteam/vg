#ifndef VG_NESTED_FLOW_CALLER_HPP_INCLUDED
#define VG_NESTED_FLOW_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
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
#include "snarl_graph.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

class SnarlGraph;
/**
 * NestedFlowCaller : DEPRECATED - Use FlowCaller with nested=true instead.
 *
 * Uses any traversals finder (ex, FlowTraversalFinder) to find
 * traversals, and calls those based on how much support they have.
 * Should work on any graph but will not report cyclic traversals.
 * This class is being replaced by FlowCaller's nested mode.
 */
class NestedFlowCaller : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
    NestedFlowCaller(const PathPositionHandleGraph& graph,
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
                     bool genotype_snarls);

    virtual ~NestedFlowCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// stuff we remember for each snarl call, to be used when genotyping its parent
    struct CallRecord {
        vector<SnarlTraversal> travs;
        vector<pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>>> genotype_by_ploidy;
        string ref_path_name;
        pair<int64_t, int64_t> ref_path_interval;
        int ref_trav_idx; // index of ref paths in CallRecord::travs
    };
    typedef map<Snarl, CallRecord, NestedCachedPackedTraversalSupportFinder::snarl_less> CallTable;

    /// update the table of calls for each child snarl (and the input snarl)
    bool call_snarl_recursive(const Snarl& managed_snarl, int ploidy,
                              const string& parent_ref_path_name, pair<size_t, size_t> parent_ref_path_interval,
                              CallTable& call_table);

    /// emit the vcf of all reference-spanning snarls
    /// The call_table needs to be completely resolved
    bool emit_snarl_recursive(const Snarl& managed_snarl, int ploidy,
                              CallTable& call_table);

    /// transform the nested allele string from something like AAC<6_10>TTT to
    /// a proper string by recursively resolving the nested snarls into alleles
    string flatten_reference_allele(const string& nested_allele, const CallTable& call_table) const;
    string flatten_alt_allele(const string& nested_allele, int allele, int ploidy, const CallTable& call_table) const;

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

    /// until we support nested snarls, cap snarl size we attempt to process
    size_t max_snarl_shallow_size = 50000;

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

    /// a hook into the snarl_caller's nested support finder
    NestedCachedPackedTraversalSupportFinder& nested_support_finder;
};

}

#endif
