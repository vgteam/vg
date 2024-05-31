#ifndef VG_DECONSTRUCTOR_HPP_INCLUDED
#define VG_DECONSTRUCTOR_HPP_INCLUDED
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <atomic>
#include "genotypekit.hpp"
#include "Variant.h"
#include "handle.hpp"
#include "traversal_finder.hpp"
#include "graph_caller.hpp"
#include "lru_cache.h"

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

// note: added VCFOutputCaller parent class from vg call bring in sorted vcf output.  it would
//       be nice to re-use more of the VCFOutputCaller code, much of which is still duplicated in
//       Deconstructor
class Deconstructor : public VCFOutputCaller {
public:

    Deconstructor();
    ~Deconstructor();

    // deconstruct the entire graph to cout.
    // Not even a little bit thread safe.
    void deconstruct(vector<string> refpaths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                     bool include_nested,
                     int context_jaccard_window,
                     bool untangle_traversals,
                     bool keep_conflicted,
                     bool strict_conflicts,
                     bool long_ref_contig,
                     double cluster_threshold = 1.0,
                     gbwt::GBWT* gbwt = nullptr,
                     bool nested_decomposition = false);
    
private:

    // initialize the vcf and get the header 
    string get_vcf_header();

    // the header needs to be initialized *before* construction for vcflib
    // but we don't know all the non-ref contigs (in nested mode) until *after*
    // construction.  end result: this hacky function to patch them in before printing
    string add_contigs_to_vcf_header(const string& vcf_header) const;
    
    // deconstruct all snarls in parallel (ie nesting relationship ignored)
    void deconstruct_graph(SnarlManager* snarl_manager);

    // deconstruct all top-level snarls in parallel
    // nested snarls are processed after their parents in the same thread
    // (same logic as vg call)
    void deconstruct_graph_top_down(SnarlManager* snarl_manager);

    // some information we pass from parent to child site when
    // doing nested deconstruction
    struct NestingInfo {
        bool has_ref;
        vector<pair<handle_t, handle_t>> child_snarls;
        PathInterval ref_path_interval;
        unordered_map<string, vector<int>> sample_to_haplotypes;
        int ref_path_allele;
    };
    
    // write a vcf record for the given site.  returns true if a record was written
    // (need to have a path going through the site)
    // the nesting_info structs are optional and used to pass reference information through nested sites...
    // the output nesting_info vector writes a record for each child snarl
    bool deconstruct_site(const handle_t& snarl_start, const handle_t& snarl_end,
                          const NestingInfo* in_nesting_info = nullptr,
                          vector<NestingInfo>* out_nesting_infos = nullptr) const;

    // get the traversals for a given site
    // this returns a combination of embedded path traversals and gbwt traversals
    // the embedded paths come first, and only they get trav_steps.
    // so you can use trav_steps.size() to find the index of the first gbwt traversal...
    void get_traversals(const handle_t& snarl_start, const handle_t& snarl_end,
                        vector<Traversal>& out_travs,
                        vector<string>& out_trav_path_names,
                        vector<pair<step_handle_t, step_handle_t>>& out_trav_steps) const;

    // this is a hack to add in * alleles -- these are haplotypes that we genotyped in the
    // parent but aren't represented in any of the traversals found in the current
    // site. *-alleles are represented as empty traversals.
    // todo: conflicts arising from alt-cycles will be able to lead to conflicting
    // results -- need to overhaul code to pass more detailed traversal information
    // from parent to child to have a chance at consistently resolving
    // star traversals are appended onto travs and trav_names
    // this funtion returns a map containing both parent and child haploty
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
    unordered_set<string> gbwt_reference_samples;
    
    // infer ploidys from gbwt when possible
    unordered_map<string, pair<int, int>> gbwt_sample_to_phase_range;

    // the ref paths
    set<string> ref_paths;

    // the off-ref paths that may be found during nested deconstruction
    // (buffered by thread)
    mutable vector<unordered_set<path_handle_t>> off_ref_paths;

    // keep track of reference samples
    set<string> ref_samples;

    // do we need to write metadata for reference contigs
    bool long_ref_contig = false;
    
    // keep track of the non-ref paths as they will be our samples
    set<string> sample_names;

    // map the path name to the sample in the vcf
    const unordered_map<string, pair<string, int>>* path_to_sample_phase;

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
    // currently implemented as handle jaccard coefficient.  So 1 means only
    // merge if identical (which is what deconstruct has always done)
    double cluster_threshold = 1.0;

    // activate the new nested decomposition mode, which is like the old include_nested
    // (which lives in vcfoutputcaller) but with more of an effort to link
    // the parent and child snarls, as well as better support for nested insertions
    bool nested_decomposition = false;
};


}
#endif
