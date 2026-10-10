#ifndef VG_CANDIDATE_FINDER_HPP_INCLUDED
#define VG_CANDIDATE_FINDER_HPP_INCLUDED

#include <limits>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

/// A set of traversals through a child snarl that are consistent with
/// a single parent allele. Multiple traversals can exist if the child
/// has internal variation within a shared region.
using TraversalSet = vector<SnarlTraversal>;

/// One TraversalSet per parent allele (index matches parent genotype).
/// For a diploid parent with genotype [0,1], element 0 contains traversals
/// consistent with parent allele 0, element 1 with parent allele 1.
using ChildTraversalSets = vector<TraversalSet>;

/**
 * Finds what a caller genotypes a site against: the reference path through the site, the site
 * oriented forward along it, and the site's candidate traversals with the reference traversal
 * among them. FlowCaller and TreeGenotyper both use it.
 */
class CandidateFinder {
public:
    /// A site ready to genotype.
    struct Site {
        /// The site, oriented forward along `ref_path_name`. A copy, not the manager's.
        Snarl snarl;
        /// The reference path the site lies on: one through both its boundaries, or else its
        /// parent's.
        string ref_path_name;
        /// Where the site lies on `ref_path_name`, as `get_ref_interval` gives it, or the parent's
        /// interval when `use_parent_interval` is set.
        tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
        /// No reference path passes through both of the site's boundaries, so it takes its
        /// parent's path and interval, and has no reference traversal of its own.
        bool use_parent_interval = false;
        /// The candidate traversals.
        vector<SnarlTraversal> travs;
        /// The index of the reference traversal in `travs`, or -1 where there is none.
        int ref_trav_idx = -1;
    };

    /// Find candidates in `graph` with `traversal_finder`, on the reference paths in
    /// `ref_path_set` (all non-alt paths where it is empty). A site is skipped unless its longest
    /// traversal is at least `allele_length_range.first` long and none is longer than
    /// `allele_length_range.second`. Nothing is owned.
    CandidateFinder(const PathPositionHandleGraph& graph,
                    TraversalFinder& traversal_finder,
                    const unordered_set<string>& ref_path_set,
                    const pair<size_t, size_t>& allele_length_range);

    /// Skip a site with more edges than this, including those of its nested sites. Zero removes
    /// the limit, which is also the default.
    void set_max_snarl_edges(size_t edges) {
        max_snarl_edges = edges ? edges : numeric_limits<size_t>::max();
    }

    /// Fill `site` for `given_snarl`, oriented as its decomposition orients it. A site on no
    /// reference path takes its parent's, `parent_ref_path_name` and `parent_ref_interval`, if it
    /// has parent traversal sets (top-down calling) or `no_reference` allows it (off-reference
    /// nested calling). Where it has parent traversal sets, the first traversal of the first
    /// non-empty set stands in for its reference traversal. Returns false, with `site` partly
    /// filled, when the site cannot be genotyped: it is one node, outside the graph, too big,
    /// outside the allele length range, or on no usable reference path.
    bool find(const Snarl& given_snarl, const string& parent_ref_path_name,
              pair<size_t, size_t> parent_ref_interval,
              const ChildTraversalSets* parent_child_trav_sets, bool no_reference,
              Site& site) const;

    /// Find all traversals through a child snarl that are consistent with a parent traversal.
    /// "Consistent" means the child's entry/exit points match what's in the parent traversal.
    /// Uses the traversal finder to enumerate all valid paths through the child.
    /// @param parent_trav The parent traversal defining entry/exit constraints
    /// @param child The child snarl to find traversals through
    /// @return Set of traversals through child, empty if parent doesn't traverse child
    TraversalSet find_child_traversal_set(const SnarlTraversal& parent_trav,
                                          const Snarl& child) const;

private:
    const PathPositionHandleGraph& graph;
    TraversalFinder& traversal_finder;
    const unordered_set<string>& ref_path_set;
    /// See set_max_snarl_edges.
    size_t max_snarl_edges = numeric_limits<size_t>::max();
    /// See the constructor.
    pair<size_t, size_t> allele_length_range;
};

}

#endif
