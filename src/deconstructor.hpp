#ifndef VG_DECONSTRUCTOR_HPP_INCLUDED
#define VG_DECONSTRUCTOR_HPP_INCLUDED
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <algorithm>
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
    void deconstruct(vector<string> refpaths, const PathPositionHandleGraph* grpah, SnarlManager* snarl_manager,
                     bool path_restricted_traversals,
                     int ploidy,
                     bool include_nested,
                     int context_jaccard_window,
                     bool untangle_traversals,
                     bool keep_conflicted,
                     bool strict_conflicts,
                     gbwt::GBWT* gbwt = nullptr);
    
private:

    // write a vcf record for the given site.  returns true if a record was written
    // (need to have a path going through the site)
    bool deconstruct_site(const Snarl* site) const;

    // convert traversals to strings.  returns mapping of traversal (offset in travs) to allele
    vector<int> get_alleles(vcflib::Variant& v,
                            const pair<vector<SnarlTraversal>,
                                       vector<pair<step_handle_t, step_handle_t>>>& path_travs,
                            int ref_path_idx,
                            const vector<bool>& use_trav,
                            char prev_char, bool use_start) const;
    
    // write traversal path names as genotypes
    void get_genotypes(vcflib::Variant& v, const vector<string>& names, const vector<int>& trav_to_allele) const;

    // given a set of traversals associated with a particular sample, select a set of size <ploidy> for the VCF
    // the highest-frequency ALT traversal is chosen
    // the bool returned is true if multiple traversals map to different alleles, more than ploidy.
    pair<vector<int>, bool> choose_traversals(const string& sample_name,
                                              const vector<int>& travs, const vector<int>& trav_to_allele,
                                              const vector<string>& trav_to_name,
                                              const vector<int>& gbwt_phases) const;

    // check to see if a snarl is too big to exhaustively traverse
    bool check_max_nodes(const Snarl* snarl) const;

    // get traversals from the exhaustive finder.  if they have nested visits, fill them in (exhaustively)
    // with node visits
    vector<SnarlTraversal> explicit_exhaustive_traversals(const Snarl* snarl) const;

    // get the path location of a given traversal out of the gbwt
    // this will be much slower than doing the same using the PathPositionGraph interface as there's no
    // underlying index. 
    tuple<bool, handle_t, size_t> get_gbwt_path_position(const SnarlTraversal& trav, const gbwt::size_type& thread) const;

    // gets a sorted node id context for a given path
    vector<nid_t> get_context(
        const pair<vector<SnarlTraversal>,
                   vector<pair<step_handle_t, step_handle_t>>>& path_travs,
        const int& trav_idx) const;

    // the underlying context-getter
    vector<nid_t> get_context(
        step_handle_t start_step,
        step_handle_t end_step) const;

    // compares node contexts
    double context_jaccard(const vector<nid_t>& target,
                           const vector<nid_t>& query) const;

    // specialization for enc_vectors
    double context_jaccard(
        const dac_vector<>& target,
        const vector<nid_t>& query) const;
    
    // toggle between exhaustive and path restricted traversal finder
    bool path_restricted = false;

    // the max ploidy we expect.
    int ploidy;

    // the graph
    const PathPositionHandleGraph* graph;

    // the snarl manager
    SnarlManager* snarl_manager;

    // the traversal finders. we always use a path traversal finder to get the reference path
    unique_ptr<PathTraversalFinder> path_trav_finder;
    // we optionally use another (exhaustive for now) traversal finder if we don't want to rely on paths
    unique_ptr<TraversalFinder> trav_finder;
    // we can also use a gbwt for traversals
    unique_ptr<GBWTTraversalFinder> gbwt_trav_finder;
    // When using the gbwt we need some precomputed information to ask about stored paths.
    unordered_set<string> gbwt_reference_samples;
    // hacky path position index for alts in the gbwt
    // we map from gbwt path id -> { map of handle -> offset } for every handle in the path
    // because child snarls are done in series, we often hit the same non-ref path consecutively
    // which makes the lru cache fairly effective
    size_t lru_size = 10; 
    vector<LRUCache<gbwt::size_type, shared_ptr<unordered_map<handle_t, size_t>>>*> gbwt_pos_caches;
    /// We need to keep track of what OMP parallelism level we made the cache
    /// list for, so we can make sure to bail out if we end up trying to use
    /// the wrong level's thread numbers.
    size_t gbwt_pos_caches_level = std::numeric_limits<size_t>::max();
    // infer ploidys from gbwt when possible
    unordered_map<string, pair<int, int>> gbwt_sample_to_phase_range;

    // the ref paths
    set<string> ref_paths;

    // keep track of reference samples
    set<string> ref_samples;
    
    // keep track of the non-ref paths as they will be our samples
    set<string> sample_names;

    // map the path name to the sample in the vcf
    const unordered_map<string, pair<string, int>>* path_to_sample_phase;

    // the sample ploidys given in the phases in our path names
    unordered_map<string, int> sample_ploidys;

    // upper limit of degree-2+ nodes for exhaustive traversal
    int max_nodes_for_exhaustive = 100;

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
};

// helpel for measuring set intersectiond and union size
template <typename T>
class count_back_inserter {
    size_t &count;
public:
    typedef void value_type;
    typedef void difference_type;
    typedef void pointer;
    typedef void reference;
    typedef std::output_iterator_tag iterator_category;
    count_back_inserter(size_t &count) : count(count) {};
    void operator=(const T &){ ++count; }
    count_back_inserter &operator *(){ return *this; }
    count_back_inserter &operator++(){ return *this; }
};

}
#endif
