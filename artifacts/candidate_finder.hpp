#ifndef VG_CANDIDATE_FINDER_HPP_INCLUDED
#define VG_CANDIDATE_FINDER_HPP_INCLUDED

#include <functional>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../src/handle.hpp"
#include "reference_path_table.hpp"
#include "../src/snarls.hpp"
#include "../src/traversal_finder.hpp"

namespace vg {

using namespace std;

/// The traversals through a child snarl that agree with one parent allele.
///
/// There can be several if the child has variation of its own.
using TraversalSet = vector<SnarlTraversal>;

/// The walks through a child site that agree with one parent allele.
///
/// There can be several if the child has variation of its own.
using WalkSet = vector<Traversal>;

/**
 * A finder of what a caller genotypes a site against.
 *
 * It finds the reference path through the site, orients the site forward along it, and finds the
 * site's candidate walks with the reference walk among them. It works on a site's bounds and
 * walks, and also offers the same for a `Snarl` and its `SnarlTraversal`s.
 */
class CandidateFinder {
public:
    /// A site ready to genotype, as its bounds and walks.
    struct WalkSite {
        /// The site's start bound, oriented into the site, forward along `ref_path_name`.
        handle_t start;
        /// The site's end bound, oriented out of the site, forward along `ref_path_name`.
        handle_t end;
        /// The reference path the site lies on: one through both its boundaries, or else its
        /// parent's.
        string ref_path_name;
        /// Where the site lies on `ref_path_name`, as `get_ref_interval` gives it, or the parent's
        /// interval when `use_parent_interval` is set.
        tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
        /// Whether the site takes its parent's path and interval, because no reference path
        /// visits both of its boundaries.
        ///
        /// Such a site has no reference walk of its own.
        bool use_parent_interval = false;
        /// The candidate walks, each from `start` to `end`.
        vector<Traversal> walks;
        /// The index of the reference walk in `walks`, or -1 where there is none.
        int ref_trav_idx = -1;
    };

    /// A site ready to genotype, as a `Snarl` and its `SnarlTraversal`s.
    struct Site {
        /// The site, oriented forward along `ref_path_name`. A copy, not the manager's.
        Snarl snarl;
        /// The reference path the site lies on: one through both its boundaries, or else its
        /// parent's.
        string ref_path_name;
        /// Where the site lies on `ref_path_name`, as `get_ref_interval` gives it, or the parent's
        /// interval when `use_parent_interval` is set.
        tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
        /// Whether the site takes its parent's path and interval, because no reference path
        /// visits both of its boundaries.
        ///
        /// Such a site has no reference traversal of its own.
        bool use_parent_interval = false;
        /// The candidate traversals.
        vector<SnarlTraversal> travs;
        /// The index of the reference traversal in `travs`, or -1 where there is none.
        int ref_trav_idx = -1;
    };

    /// Find candidates in `graph` with `traversal_finder`, on the reference paths in
    /// `ref_path_set`, or on all non-alt paths if it is empty.
    ///
    /// Skip a site unless its longest traversal is at least `allele_length_range.first` long and
    /// none is longer than `allele_length_range.second`. Nothing is owned.
    CandidateFinder(const PathPositionHandleGraph& graph,
                    SiteTraversalFinder& traversal_finder,
                    const unordered_set<string>& ref_path_set,
                    const pair<size_t, size_t>& allele_length_range);

    /// Find candidates in `graph` with `traversal_finder`, on the reference paths in
    /// `ref_paths`, or on all non-alt paths if it is empty.
    ///
    /// Skip a site unless its longest traversal is at least `allele_length_range.first` long and
    /// none is longer than `allele_length_range.second`. Nothing is owned.
    CandidateFinder(const PathPositionHandleGraph& graph,
                    SiteTraversalFinder& traversal_finder,
                    const multipass::ReferencePathTable& ref_paths,
                    const pair<size_t, size_t>& allele_length_range);

    /// Skip a site with more edges than this, including those of its nested sites.
    ///
    /// Zero removes the limit, which is also the default.
    void set_max_snarl_edges(size_t edges);

    /**
     * Fill `site` for the site between `start`, oriented into it, and `end`, oriented out of it.
     *
     * A site on no reference path takes its parent's, `parent_ref_path_name` and
     * `parent_ref_interval`, if it has parent walk sets (top-down calling) or `no_reference`
     * allows it (off-reference nested calling). With parent walk sets, the first walk of the first
     * nonempty set stands in for its reference walk.
     *
     * @param parent_child_walk_sets If not null, one WalkSet per allele of the parent's genotype,
     * in the order of the genotype.
     * @return False, with `site` partly filled, if the site cannot be genotyped. That is when it
     * is one node, outside the graph, too big, outside the allele length range, or on no usable
     * reference path.
     */
    bool find(const handle_t& start, const handle_t& end, const string& parent_ref_path_name,
              pair<size_t, size_t> parent_ref_interval,
              const vector<WalkSet>* parent_child_walk_sets, bool no_reference,
              WalkSite& site) const;

    /// Fill `site` for `given_snarl`, as `find` does for its bounds.
    bool find(const Snarl& given_snarl, const string& parent_ref_path_name,
              pair<size_t, size_t> parent_ref_interval,
              const vector<TraversalSet>* parent_child_trav_sets, bool no_reference,
              Site& site) const;

    /// Find the walks through the child site between `child_start` and `child_end` that agree
    /// with the parent walk `parent_walk`.
    ///
    /// @return Every walk the traversal finder finds through the child, if `parent_walk` visits
    /// both of the child's boundary nodes, or else none.
    WalkSet find_child_walk_set(const Traversal& parent_walk, const handle_t& child_start,
                                const handle_t& child_end) const;

    /// Find the traversals through `child` that agree with `parent_trav`, as
    /// `find_child_walk_set` does.
    TraversalSet find_child_traversal_set(const SnarlTraversal& parent_trav,
                                          const Snarl& child) const;

private:
    /// Return true if the site between `start` and `end` has more edges than the cap.
    bool too_big(handle_t start, handle_t end) const;

    /// Get the names of the usable reference paths that visit `handle`, sorted.
    set<string> reference_paths_on(handle_t handle) const;

    /// Get the names of the usable reference paths that visit both `start` and `end`, sorted.
    vector<string> common_reference_paths(handle_t start, handle_t end) const;

    /// Trace the reference path through the site from `start` to `end` over `ref_interval` into
    /// `ref_walk`.
    ///
    /// @return False, after a warning, if the path ends before the end of the site.
    bool reference_walk(const handle_t& start, const handle_t& end,
                        const tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t>&
                            ref_interval,
                        Traversal& ref_walk) const;

    /// Return true if the longest of `walks` is within the allele length range.
    bool in_allele_length_range(const vector<Traversal>& walks) const;

    const PathPositionHandleGraph& graph;
    SiteTraversalFinder& traversal_finder;
    /// Tests whether a non-alt path, by name, is one a site can be placed on.
    function<bool(const string&)> usable_path;
    /// The most edges a site may have, including those of its nested sites, or the largest
    /// `size_t` for no limit.
    size_t max_snarl_edges = numeric_limits<size_t>::max();
    /// The bounds on the length, in bases, of the longest traversal of a site.
    pair<size_t, size_t> allele_length_range;
};

}

#endif
