#include "handle.hpp"
#include "traversal_finder.hpp"

#pragma once


/** \file
* Utilities for finding and clustering similar snarl traversals
*/
namespace vg {
using namespace std;

// 
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

// compute the jaccard coefficient (|intersection| / |union|) of two sets
// target and query should be sorted vectors (or something equivalent)
template <typename T, typename U>
inline double jaccard_coefficient(const T& target, const U& query) {
    size_t isec_size = 0;
    std::set_intersection(target.begin(), target.end(),
                          query.begin(), query.end(),
                          count_back_inserter<typename T::value_type>(isec_size));
    size_t union_size = 0;
    std::set_union(target.begin(), target.end(),
                   query.begin(), query.end(),
                   count_back_inserter<typename T::value_type>(union_size));
    return (double)isec_size / (double)union_size;
}

// specialized version that weights jaccard by node lenths. 
double weighted_jaccard_coefficient(const PathHandleGraph* graph, const multiset<handle_t>& target, const multiset<handle_t>& query);

// the information needed from the parent traversal in order to
// genotype a child traversal
struct ParentGenotypeInfo {
    Traversal ref_traversal;
    pair<step_handle_t, step_handle_t> ref_steps;
    unordered_map<string, int64_t> sample_to_ploidy; // to check/enforce consistency
};


/// sort the traversals, putting the reference first then using names
/// traversals masked out by use_traversal will be filrtered out entirely
/// (so the output vector may be smaller than the input...)
vector<int> get_traversal_order(const PathHandleGraph* graph,
                                const vector<Traversal>& traversals,
                                const vector<string>& trav_path_names,
                                const vector<int>& ref_travs,
                                const vector<bool>& use_traversal);


/// cluster the traversals. The algorithm is:
/// - visit traversals in provided order
///   - if the traversal is <= min_jaccard away from the reference traversal of cluster, add to cluster
///   - else start a new cluster, with the given traversal as a reference
/// note that traversal_order can specify a subset of traversals
/// out_info are Simlarity/length-delta pairs comparing the traversal to its cluster reference
/// (if the traversal wasn't in traversal_order, it'll get -1)
/// if child snarls are given, then we make sure that every child snarl is contained in a reference
/// traversal (this guarantees the snarl nesting structure is exactly represented in the vcf alleles!)
/// if child snarls are given, then we also fill in the mapping of each child snarl to its
/// "reference" traversals -- ie first traversal in a cluster that contains both its endpoints
vector<vector<int>> cluster_traversals(const PathHandleGraph* graph,
                                       const vector<Traversal>& traversals,
                                       const vector<int>& traversal_order,
                                       const vector<pair<handle_t, handle_t>>& child_snarls,
                                       double min_jaccard,
                                       vector<pair<double, int64_t>>& out_info,
                                       vector<int>& out_child_snarl_to_trav);

/// assign a list of child snarls to traversals that fully conrtain them
/// the output is a list (for each traversal) of each child snarl that's contained in it
/// so otuput[i] = {x,y,z} means that child_snarls[x],[y],[z] are in traversals[i]
vector<vector<int>> assign_child_snarls_to_traversals(const PathHandleGraph* graph,
                                                      const vector<Traversal>& traversals,
                                                      const vector<pair<handle_t, handle_t>>& child_snarls);

}
