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

namespace vg {

using namespace std;

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager,
                ostream& out_stream = cout);

    virtual ~GraphCaller();

    /// Run call_snarl() on every top-level snarl in the manager.
    /// For any that return false, try the children, etc. (when recurse_on_fail true)
    /// Snarls are processed in parallel
    virtual void call_top_level_snarls(bool recurse_on_fail = true);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;
   
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// Where all output written
    ostream& out_stream;
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
    void write_variants(ostream& out_stream) const;
    
protected:

    /// get the interval of a snarl from our reference path using the PathPositionHandleGraph interface
    /// the bool is true if the snarl's backward on the path
    tuple<size_t, size_t, bool, step_handle_t, step_handle_t> get_ref_interval(const PathPositionHandleGraph& graph, const Snarl& snarl,
                                                                               const string& ref_path_name) const;
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;

    /// output buffers (1/thread) (for sorting)
    mutable vector<vector<vcflib::Variant>> output_variants;
};
    
/**
 * VCFGenotyper : Genotype variants in a given VCF file
 */
class VCFGenotyper : public GraphCaller, public VCFOutputCaller {
public:
    VCFGenotyper(const PathHandleGraph& graph,
                 SnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 vcflib::VariantCallFile& variant_file,
                 const string& sample_name,
                 const vector<string>& ref_paths = {},
                 FastaReference* ref_fasta = nullptr,
                 FastaReference* ins_fasta = nullptr,
                 ostream& out_stream = cout);

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
                 const vector<size_t>& ref_path_offsets = {});

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

    /// print a vcf variant 
    void emit_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const unique_ptr<SnarlCaller::CallInfo>& call_info, const string& ref_path_name) const;

    /// check if a site can be handled by the RepresentativeTraversalFinder
    bool is_traversable(const Snarl& snarl);

    /// look up a path index for a site and return its name too
    pair<string, PathIndex*> find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const;

    /// clean up the alleles to not share common prefixes / suffixes
    void flatten_common_allele_ends(vcflib::Variant& variant, bool backward) const;

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
 * FlowCaller : Uses the FlowTraversal finder to find best-supported
 * traversals, and calls those.  should work on any graph but will not
 * report cyclic traversals.  Does not (yet, anyway) support nested
 * calling, so the entire site is processes in one shot. 
 * Designed to replace LegacyCaller, as it should miss fewer obviously
 * good traversals, and is not dependent on old protobuf-based structures. 
 */
class FlowCaller : public GraphCaller, public VCFOutputCaller {
public:
    FlowCaller(const PathPositionHandleGraph& graph,
               SupportBasedSnarlCaller& snarl_caller,
               SnarlManager& snarl_manager,
               const string& sample_name,
               size_t max_traverals = 100,
               const vector<string>& ref_paths = {},
               const vector<size_t>& ref_path_offsets = {},
               ostream& out_stream = cout);
   
    virtual ~FlowCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    // TODO:
    // these methods can and should be merged with legacy caller, maybe by pushing up to VCFOutputCaller

    /// print a vcf variant 
    void emit_variant(const Snarl& snarl, int ref_trav_idx, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const unique_ptr<SnarlCaller::CallInfo>& call_info, const string& ref_path_name) const;

    /// clean up the alleles to not share common prefixes / suffixes
    void flatten_common_allele_ends(vcflib::Variant& variant, bool backward) const;

protected:

    /// the graph
    const PathPositionHandleGraph& graph;

    /// the traversal finder
    FlowTraversalFinder* traversal_finder;

    /// keep track of the reference paths
    vector<string> ref_paths;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;

};



}

#endif
