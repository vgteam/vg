#ifndef VG_GRAPH_CALLER_HPP_INCLUDED
#define VG_GRAPH_CALLER_HPP_INCLUDED

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

namespace vg {

using namespace std;

using vg::io::AlignmentEmitter;

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:

    enum RecurseType { RecurseOnFail, RecurseAlways, RecurseNever };
   
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager);

    virtual ~GraphCaller();

    /// Run call_snarl() on every top-level snarl in the manager.
    /// For any that return false, try the children, etc. (when recurse_on_fail true)
    /// Snarls are processed in parallel
    virtual void call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type = RecurseOnFail);

    /// For every chain, cut it up into pieces using max_edges and max_trivial to cap the size of each piece
    /// then make a fake snarl for each chain piece and call it.  If a fake snarl fails to call,
    /// It's child chains will be recursed on (if selected)_
    virtual void call_top_level_chains(const HandleGraph& graph,
                                       size_t max_edges,
                                       size_t max_trivial,
                                       RecurseType recurise_type = RecurseOnFail);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

    /// toggle progress messages
    void set_show_progress(bool show_progress);

protected:

    /// Break up a chain into bits that we want to call using size heuristics
    vector<Chain> break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial);
    
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// Toggle progress messages
    bool show_progress;
};

/**
 * Helper class that vcf writers can inherit from to for some common code to output sorted VCF
 */
class VCFOutputCaller {
public:
    VCFOutputCaller(const string& sample_name);

    virtual ~VCFOutputCaller();

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const;

    /// Add a variant to our buffer
    void add_variant(vcflib::Variant& var) const;

    /// Sort then write variants in the buffer
    /// snarl_manager needed if include_nested is true
    void write_variants(ostream& out_stream, const SnarlManager* snarl_manager = nullptr);

    /// Run vcffixup from vcflib
    void vcf_fixup(vcflib::Variant& var) const;

    /// Add a translation map
    void set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation);

    /// Assume writing nested snarls is enabled
    void set_nested(bool nested);
    
protected:

    /// add a traversal to the VCF info field in the format of a GFA W-line or GAF path
    void add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele,
                                 const Traversal& trav, bool reversed, bool one_based) const;
    /// legacy version of above
    void add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav, bool reversed, bool one_based) const;
    
    
    /// convert a traversal into an allele string
    string trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const;
    
    /// print a vcf variant
    void emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                      function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string = nullptr);

    /// get the interval of a snarl from our reference path using the PathPositionHandleGraph interface
    /// the bool is true if the snarl's backward on the path
    /// first returned value -1 if no traversal found 
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(const PathPositionHandleGraph& graph, const Snarl& snarl,
                                                                                 const string& ref_path_name) const;

    /// used for making gaf traversal names
    pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                           int64_t ref_path_offset) const;

    /// clean up the alleles to not share common prefixes / suffixes
    /// if len_override given, just do that many bases without thinking
    void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const;

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;

    /// do the opposite of above
    /// So a string that looks like AACT(>12<17)TTT would invoke the callback three times with
    /// ("AACT", Snarl), ("", Snarl(12,-17)), ("TTT", Snarl(12,-17))
    /// The parameters are to be treated as unions:  A sequence fragment if non-empty, otherwise a snarl
    void scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const;

    // update the PS and LV tags in the output buffer (called in write_variants if include_nested is true)
    void update_nesting_info_tags(const SnarlManager* snarl_manager);
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;

    /// output buffers (1/thread) (for sorting)
    /// variants stored as strings (and position key pairs) because vcflib::Variant in-memory struct so huge
    mutable vector<vector<pair<pair<string, size_t>, string>>> output_variants;

    /// print up to this many uncalled alleles when doing ref-genotpes in -a mode
    size_t max_uncalled_alleles = 5;

    // optional node translation to apply to snarl names in variant IDs
    const unordered_map<nid_t, pair<string, size_t>>* translation;

    // need to write LV/PS info tags
    bool include_nested;
};

/**
 * Helper class for outputing snarl traversals as GAF
 */
class GAFOutputCaller {
public:
    /// The emitter object is created and owned by external forces
    GAFOutputCaller(AlignmentEmitter* emitter, const string& sample_name, const vector<string>& ref_paths,
                    size_t trav_padding);
    virtual ~GAFOutputCaller();

    /// print the GAF traversals
    void emit_gaf_traversals(const PathHandleGraph& graph, const string& snarl_name,
                             const vector<SnarlTraversal>& travs,
                             int64_t ref_trav_idx,
                             const string& ref_path_name, int64_t ref_path_position,
                             const TraversalSupportFinder* support_finder = nullptr);

    /// print the GAF genotype
    void emit_gaf_variant(const PathHandleGraph& graph, const string& snarl_name,
                          const vector<SnarlTraversal>& travs,
                          const vector<int>& genotype,
                          int64_t ref_trav_idx,
                          const string& ref_path_name, int64_t ref_path_position,
                          const TraversalSupportFinder* support_finder = nullptr);
    
    /// pad a traversal with (first found) reference path, adding up to trav_padding to each side
    SnarlTraversal pad_traversal(const PathHandleGraph& graph, const SnarlTraversal& trav) const;
    
protected:
    
    AlignmentEmitter* emitter;

    /// Sample name
    string gaf_sample_name;

    /// Add padding from reference paths to traversals to make them at least this long
    /// (only in emit_gaf_traversals(), not emit_gaf_variant)
    size_t trav_padding = 0;

    /// Reference paths are used to pad out traversals.  If there are none, then first path found is used
    unordered_set<string> ref_paths;

};

/**
 * VCFGenotyper : Genotype variants in a given VCF file
 */
class VCFGenotyper : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
    VCFGenotyper(const PathHandleGraph& graph,
                 SnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 vcflib::VariantCallFile& variant_file,
                 const string& sample_name,
                 const vector<string>& ref_paths,
                 const vector<int>& ref_path_ploidies,
                 FastaReference* ref_fasta,
                 FastaReference* ins_fasta,
                 AlignmentEmitter* aln_emitter,
                 bool traversals_only,
                 bool gaf_output,
                 size_t trav_padding);

    virtual ~VCFGenotyper();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// get path positions bounding a set of variants
    tuple<string, size_t, size_t>  get_ref_positions(const vector<vcflib::Variant*>& variants) const;

    /// munge out the contig lengths from the VCF header
    virtual unordered_map<string, size_t> scan_contig_lengths() const;

protected:

    /// the graph
    const PathHandleGraph& graph;

    /// input VCF to genotype, must have been loaded etc elsewhere
    vcflib::VariantCallFile& input_vcf;

    /// traversal finder uses alt paths to map VCF alleles from input_vcf
    /// back to traversals in the snarl
    VCFTraversalFinder traversal_finder;

    /// toggle whether to genotype or just output the traversals
    bool traversals_only;

    /// toggle whether to output vcf or gaf
    bool gaf_output;

    /// the ploidies
    unordered_map<string, int> path_to_ploidy;
};


/**
 * LegacyCaller : Preserves (most of) the old vg call logic by using 
 * the RepresentativeTraversalFinder to recursively find traversals
 * through arbitrary sites.   
 */
class LegacyCaller : public GraphCaller, public VCFOutputCaller {
public:
    LegacyCaller(const PathPositionHandleGraph& graph,
                 SupportBasedSnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 const string& sample_name,
                 const vector<string>& ref_paths = {},
                 const vector<size_t>& ref_path_offsets = {},
                 const vector<int>& ref_path_ploidies = {});

    virtual ~LegacyCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// recursively genotype a snarl
    /// todo: can this be pushed to a more generic class? 
    pair<vector<SnarlTraversal>, vector<int>> top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy,
                                                                const string& ref_path_name, pair<size_t, size_t> ref_interval) const;
    
    /// we need the reference traversal for VCF, but if the ref is not called, the above method won't find it. 
    SnarlTraversal get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const;

    /// re-genotype output of top_down_genotype.  it may give slightly different results as
    /// it's working with fully-defined traversals and can exactly determine lengths and supports
    /// it will also make sure the reference traversal is in the beginning of the output
    tuple<vector<SnarlTraversal>, vector<int>, unique_ptr<SnarlCaller::CallInfo>> re_genotype(const Snarl& snarl,
                                                                                              TraversalFinder& trav_finder,
                                                                                              const vector<SnarlTraversal>& in_traversals,
                                                                                              const vector<int>& in_genotype,
                                                                                              int ploidy,
                                                                                              const string& ref_path_name,
                                                                                              pair<size_t, size_t> ref_interval) const;

    /// check if a site can be handled by the RepresentativeTraversalFinder
    bool is_traversable(const Snarl& snarl);

    /// look up a path index for a site and return its name too
    pair<string, PathIndex*> find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const;

protected:

    /// the graph
    const PathPositionHandleGraph& graph;
    /// non-vg inputs are converted into vg as-needed, at least until we get the
    /// traversal finding ported
    bool is_vg;

    /// The old vg call traversal finder.  It is fairly efficient but daunting to maintain.
    /// We keep it around until a better replacement is implemented.  It is *not* compatible
    /// with the Handle Graph API because it relise on PathIndex.  We convert to VG as
    /// needed in order to use it. 
    RepresentativeTraversalFinder* traversal_finder;
    /// Needed by above (only used when working on vg inputs -- generated on the fly otherwise)
    vector<PathIndex*> path_indexes;

    /// keep track of the reference paths
    vector<string> ref_paths;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;

    /// keep track of ploidies in the reference paths
    map<string, int> ref_ploidies;

    /// Tuning

    /// How many nodes should we be willing to look at on our path back to the
    /// primary path? Keep in mind we need to look at all valid paths (and all
    /// combinations thereof) until we find a valid pair.
    int max_search_depth = 1000;
    /// How many search states should we allow on the DFS stack when searching
    /// for traversals?
    int max_search_width = 1000;
    /// What's the maximum number of bubble path combinations we can explore
    /// while finding one with maximum support?
    size_t max_bubble_paths = 100;

};


/**
 * FlowCaller : Uses any traversals finder (ex, FlowTraversalFinder) to find 
 * traversals, and calls those based on how much support they have.  
 * Should work on any graph but will not
 * report cyclic traversals.  Does not (yet, anyway) support nested
 * calling, so the entire site is processes in one shot. 
 * Designed to replace LegacyCaller, as it should miss fewer obviously
 * good traversals, and is not dependent on old protobuf-based structures. 
 */
class FlowCaller : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
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
   
    virtual ~FlowCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

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

    /// until we support nested snarls, cap snarl size we attempt to process
    size_t max_snarl_edges = 10000;

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
};

class SnarlGraph;

/**
 * FlowCaller : Uses any traversals finder (ex, FlowTraversalFinder) to find 
 * traversals, and calls those based on how much support they have.  
 * Should work on any graph but will not
 * report cyclic traversals.  
 *
 * todo: this is a generalization of FlowCaller and should be able to replace it entirely after testing
 *       to get rid of duplicated code. 
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


/** Simplification of a NetGraph that ignores chains.  It is designed only for
    traversal finding.  Todo: generalize NestedFlowCaller to the point where we 
    can remove this and use NetGraph instead */
class SnarlGraph : virtual public HandleGraph {
public:
    // note: can only deal with one snarl "level" at a time
    SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls);

    // go from node to snarl (first val false if not a snarl)
    pair<bool, handle_t> node_to_snarl(handle_t handle) const;

    // go from edge to snarl (first val false if not a virtual edge)
    tuple<bool, handle_t, edge_t> edge_to_snarl_edge(edge_t edge) const;

    // replace a snarl node with an actual snarl in the traversal
    void embed_snarl(Visit& visit);
    void embed_snarls(SnarlTraversal& traversal);

    // replace a refpath through the snarl with the actual snarl in the traversal
    // todo: this is a bed of a hack
    void embed_ref_path_snarls(SnarlTraversal& traversal);

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface (which is all identical to backing graph)
    ////////////////////////////////////////////////////////////////////////////
    bool has_node(nid_t node_id) const;
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    nid_t get_id(const handle_t& handle) const;
    bool get_is_reverse(const handle_t& handle) const;
    handle_t flip(const handle_t& handle) const;
    size_t get_length(const handle_t& handle) const;
    std::string get_sequence(const handle_t& handle) const;    
    size_t get_node_count() const;
    nid_t min_node_id() const;
    nid_t max_node_id() const;
    
protected:

    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// this is the only function that's changed to do anything different from the backing graph:
    /// it is changed to "pass through" snarls by pretending there are edges from into snarl starts out of ends and
    /// vice versa.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;    

    /// the backing graph
    const HandleGraph* backing_graph;

    /// the snarl manager
    SnarlManager& snarl_manager;

    /// the snarls (indexed both ways).  flag is true for original orientation
    unordered_map<handle_t, pair<handle_t, bool>> snarls;
};


}

#endif
