#ifndef VG_DECONSTRUCTOR_HPP_INCLUDED
#define VG_DECONSTRUCTOR_HPP_INCLUDED
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <atomic>
#include <memory>
#include "genotypekit.hpp"
#include "Variant.h"
#include "handle.hpp"
#include "traversal_finder.hpp"
#include "graph_caller.hpp"
#include "sorted_vcf_writer.hpp"
#include "lru_cache.h"

#include <gbwtgraph/utils.h>

/** \file
* Deconstruct is getting rewritten.
* New functionality:
* -Detect superbubbles and bubbles
* -Fix command line interface.
* -harmonize on XG / raw graph (i.e. deprecate index)
* -Use unroll/DAGify if needed to avoid cycles

** Much of this is taken from Brankovic's
** "Linear-Time Superbubble Identification Algorithm for Genome Assembly"
*/
namespace vg{
using namespace std;

// note: writes through a SortedVCFWriter, as vg call does, for sorted vcf output.  it would
//       be nice to re-use more of vg call's VCF code, much of which is still duplicated in
//       Deconstructor
class Deconstructor {
public:

    Deconstructor();
    ~Deconstructor();

    /// The sorted VCF the records go to.
    SortedVCFWriter vcf_writer;

    // deconstruct the entire graph to cout.
    // Not even a little bit thread safe.
    //
    // With `include_nested`, `nesting_sites` must be the snarls of `snarl_manager` as a snarl
    // decomposition, which the nesting INFO tags are read from.
    void deconstruct(vector<string> refpaths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                     const SnarlDecomposition* nesting_sites,
                     bool include_nested,
                     int context_jaccard_window,
                     bool untangle_traversals,
                     bool keep_conflicted,
                     bool strict_conflicts,
                     bool long_ref_contig,
                     double cluster_threshold = 1.0,
                     int64_t cluster_min_allele_len = 0,
                     gbwt::GBWT* gbwt = nullptr,
                     bool star_allele = false);
    
private:

    // Find the samples, from the paths that are not references, and the ploidy of each sample.
    void find_samples();

    // Declare in the header of `vcf_writer` the fields that the records have, other than GT and
    // the nesting tags. The samples must have been found.
    void declare_fields();

    // Return the reference contigs, each named as the CHROM column of its records names it. The
    // length of a contig is the furthest that any of its reference paths reaches, counting the
    // start of the path's subrange.
    vector<SortedVCFWriter::ReferenceContig> reference_contigs() const;

    // Test if a path is the reference we're deconstructing against, seen through the
    // other half of the base/gref pair (a gref path when deconstructing against the
    // base reference, or the base path of a selected gref reference).  Such a path is
    // not a sample: it would genotype as all-reference and inflate AC/AF/AN/NS.
    bool is_other_reference_view(const string& path_name) const;

    // Samples that exist only as the other view (base or gref) of a selected reference.
    set<string> other_ref_samples;

    
    // Simple flat processing: collect all snarls then process in parallel.
    // More memory-efficient when star alleles aren't needed.
    void deconstruct_graph(SnarlManager* snarl_manager);

    // Top-down processing: process snarls level by level, passing context to children.
    // Required for star allele detection (context tracks parent haplotypes).
    void deconstruct_graph_top_down(SnarlManager* snarl_manager);

    // Context passed from parent to child site for nested deconstruction.
    // Used for star allele detection: samples in parent but not in child get * allele.
    // Note: context is only populated when star_allele is enabled.
    // Note: context is NOT passed when parent has multiple reference traversals (cycles).
    struct ChildContext {
        // Map from sample name to vector of haplotype phases that traversed the parent snarl.
        // Child compares against its own traversals to detect star alleles.
        // Format matches what add_star_traversals() expects.
        // Uses shared_ptr so all children of a snarl share the same map (memory optimization).
        std::shared_ptr<const unordered_map<string, vector<int>>> sample_to_haplotypes;

        // Check if context is populated
        bool empty() const { return !sample_to_haplotypes || sample_to_haplotypes->empty(); }
    };

    // write a vcf record for the given site.  returns true if a record was written
    // (need to have a path going through the site)
    // in_context: optional context from parent (for star allele detection)
    // out_child_contexts: if provided, populated with contexts for each child snarl
    //                     (only populated when site has exactly one reference traversal)
    bool deconstruct_site(const handle_t& snarl_start, const handle_t& snarl_end,
                          const ChildContext* in_context = nullptr,
                          vector<ChildContext>* out_child_contexts = nullptr);

    // get the traversals for a given site
    // this returns a combination of embedded path traversals and gbwt traversals
    // the embedded paths come first, and only they get trav_steps.
    // so you can use trav_steps.size() to find the index of the first gbwt traversal...
    void get_traversals(const handle_t& snarl_start, const handle_t& snarl_end,
                        vector<Traversal>& out_travs,
                        vector<string>& out_trav_path_names,
                        vector<pair<step_handle_t, step_handle_t>>& out_trav_steps) const;

    // Add star (*) alleles for haplotypes that were genotyped in the parent but don't
    // have a traversal through the current site. Star alleles are represented as empty
    // traversals. Returns a map of sample name to haplotype phases (both parent and child).
    unordered_map<string, vector<int>> add_star_traversals(vector<Traversal>& travs,
                                                           vector<string>& trav_names,
                                                           vector<vector<int>>& trav_clusters,
                                                           vector<pair<double, int64_t>>& trav_cluster_info,
                                                           const unordered_map<string, vector<int>>& parent_haplotypes) const;

    // convert traversals to strings.  returns mapping of traversal (offset in travs) to allele
    vector<int> get_alleles(vcflib::Variant& v,
                            const vector<Traversal>& travs,
                            const vector<pair<step_handle_t, step_handle_t>>& trav_steps,
                            int ref_path_idx,
                            const vector<vector<int>>& trav_clusters,
                            char prev_char, bool use_start) const;
    
    // write traversal path names as genotypes
    void get_genotypes(vcflib::Variant& v, const vector<string>& names, const vector<int>& trav_to_allele,
                       const vector<pair<double, int64_t>>& trav_to_cluster_info) const;

    // given a set of traversals associated with a particular sample, select a set of size <ploidy> for the VCF
    // the highest-frequency ALT traversal is chosen
    // the bool returned is true if multiple traversals map to different alleles, more than ploidy.
    pair<vector<int>, bool> choose_traversals(const string& sample_name,
                                              const vector<int>& travs, const vector<int>& trav_to_allele,
                                              const vector<string>& trav_to_name,
                                              const vector<int>& gbwt_phases) const;

    // the underlying context-getter
    vector<nid_t> get_context(
        step_handle_t start_step,
        step_handle_t end_step) const;
    
    // the graph
    const PathPositionHandleGraph* graph;

    // the gbwt
    gbwt::GBWT* gbwt;

    // the traversal finders. we always use a path traversal finder to get the reference path
    unique_ptr<PathTraversalFinder> path_trav_finder;
    // we can also use a gbwt for traversals
    unique_ptr<GBWTTraversalFinder> gbwt_trav_finder;
    // When using the gbwt we need some precomputed information to ask about stored paths.
    gbwtgraph::sample_name_set gbwt_reference_samples;
    
    // infer ploidys from gbwt when possible
    unordered_map<string, pair<int, int>> gbwt_sample_to_phase_range;

    // the ref paths
    set<string> ref_paths;

    // keep track of reference samples
    set<string> ref_samples;

    // do we need to write metadata for reference contigs
    bool long_ref_contig = false;
    
    // keep track of the non-ref paths as they will be our samples
    set<string> sample_names;

    // the sample ploidys given in the phases in our path names
    unordered_map<string, int> sample_ploidys;

    // target window size for determining the correct reference position for allele traversals with path jaccard
    int path_jaccard_window = 10000;

    // should we add positional untangling of traversals in the AP field
    bool untangle_allele_traversals = false;

    // should we be strict about flagging and removing conflicted phases?
    bool strict_conflict_checking = false;

    // show path info mapping paths to genotypes (very verbose)
    bool show_path_info = false;

    // should we keep conflicted genotypes or not
    bool keep_conflicted_genotypes = false;

    // used to merge together similar traversals (to keep allele counts down)
    // currently a length-weighted similarity (see weighted_traversal_similarity).  So 1 means only
    // merge if identical (which is what deconstruct has always done)
    double cluster_threshold = 1.0;

    // When cluster_threshold is below 1, apply it only at sites whose core length is at least this
    // many bp.  At 0, cluster_threshold applies at every site.  The core length is the length of
    // the longest allele once the prefix and suffix that every allele other than "*" shares are
    // stripped.
    int64_t cluster_min_allele_len = 0;

    // use *-alleles to represent haplotypes that span the parent site but don't
    // traverse the current nested site (e.g., a deletion that spans a nested SNP)
    // only works with include_nested
    bool star_allele = false;

    // deconstruct the nested snarls too, and write the nesting INFO tags
    bool include_nested = false;
};


}
#endif
