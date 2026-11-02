#include <algorithm>
#include <limits>

#include "candidate_finder.hpp"
#include "../src/graph_caller.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg {

void CandidateFinder::set_max_snarl_edges(size_t edges) {
    max_snarl_edges = edges ? edges : numeric_limits<size_t>::max();
}

CandidateFinder::CandidateFinder(const PathPositionHandleGraph& graph,
                                 SiteTraversalFinder& traversal_finder,
                                 const unordered_set<string>& ref_path_set,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph), traversal_finder(traversal_finder),
    usable_path([&ref_path_set](const string& name) {
        return ref_path_set.empty() || ref_path_set.contains(name);
    }),
    allele_length_range(allele_length_range) {
}

CandidateFinder::CandidateFinder(const PathPositionHandleGraph& graph,
                                 SiteTraversalFinder& traversal_finder,
                                 const multipass::ReferencePathTable& ref_paths,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph), traversal_finder(traversal_finder),
    usable_path([&ref_paths](const string& name) {
        return ref_paths.empty() || ref_paths.contains(name);
    }),
    allele_length_range(allele_length_range) {
}

/// Convert a SnarlTraversal of node visits to a walk.
static Traversal walk_of_traversal(const HandleGraph& graph, const SnarlTraversal& trav) {
    Traversal walk;
    walk.reserve(trav.visit_size());
    for (const Visit& visit : trav.visit()) {
        walk.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return walk;
}

/// Convert a walk to a SnarlTraversal of node visits.
static SnarlTraversal traversal_of_walk(const HandleGraph& graph, const Traversal& walk) {
    SnarlTraversal trav;
    for (const handle_t& handle : walk) {
        Visit* visit = trav.add_visit();
        visit->set_node_id(graph.get_id(handle));
        visit->set_backward(graph.get_is_reverse(handle));
    }
    return trav;
}

/// Describe the site between `start` and `end` for a warning.
static string describe_site(const HandleGraph& graph, const handle_t& start, const handle_t& end) {
    return "from node " + std::to_string(graph.get_id(start))
           + (graph.get_is_reverse(start) ? " reverse" : " forward") + " to node "
           + std::to_string(graph.get_id(end))
           + (graph.get_is_reverse(end) ? " reverse" : " forward");
}

bool CandidateFinder::too_big(handle_t start, handle_t end) const {
    // Counting the edges walks the whole site, so skip it when there is no cap.
    return max_snarl_edges != numeric_limits<size_t>::max()
           && site_contents(graph, start, end, false).second.size() > max_snarl_edges;
}

set<string> CandidateFinder::reference_paths_on(handle_t handle) const {
    set<string> names;
    graph.for_each_step_on_handle(handle, [&](step_handle_t step_handle) {
        string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
        if (!Paths::is_alt(name) && usable_path(name)) {
            names.insert(name);
        }
        return true;
    });
    return names;
}

vector<string> CandidateFinder::common_reference_paths(handle_t start, handle_t end) const {
    const set<string> start_path_names = reference_paths_on(start);
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        end_path_names = reference_paths_on(end);
    }
    // Take the full intersection, rather than the first common path, so that the choice from it
    // does not depend on the order of iteration, which could change between implementations or
    // runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));
    return common_names;
}

/// Choose the reference path to call a site on from the sorted, nonempty `common_names`.
///
/// Prefer the first path outside the gref namespace, or else the first path. A gref copy of the
/// reference sorts before the path it was copied from (gref_x < x), so test the whole gref
/// namespace, not just the fragment suffix.
static string preferred_reference(const vector<string>& common_names) {
    for (const string& name : common_names) {
        if (!GrefCover::is_gref_derived(name)) {
            return name;
        }
    }
    return common_names.front();
}

bool CandidateFinder::reference_walk(const handle_t& start, const handle_t& end,
                                     const tuple<int64_t, int64_t, bool, step_handle_t,
                                                 step_handle_t>& ref_interval,
                                     Traversal& ref_walk) const {
    step_handle_t cur_step = get<3>(ref_interval);
    step_handle_t last_step = get<4>(ref_interval);
    if (get<2>(ref_interval)) {
        std::swap(cur_step, last_step);
    }
    const bool start_backwards =
        graph.get_is_reverse(start) != graph.get_is_reverse(graph.get_handle_of_step(cur_step));

    while (true) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        ref_walk.push_back(start_backwards ? graph.flip(cur_handle) : cur_handle);
        if (graph.get_id(cur_handle) == graph.get_id(end)) {
            break;
        }
        const bool path_backward = get<2>(ref_interval);
        if (path_backward ? !graph.has_previous_step(cur_step) : !graph.has_next_step(cur_step)) {
            cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace "
                    "reference path through the site "
                 << describe_site(graph, start, end) << endl;
            return false;
        }
        cur_step = path_backward ? graph.get_previous_step(cur_step) : graph.get_next_step(cur_step);
    }
    assert(ref_walk.front() == start && ref_walk.back() == end);
    return true;
}

bool CandidateFinder::in_allele_length_range(const vector<Traversal>& walks) const {
    if (allele_length_range.first == 0
        && allele_length_range.second == numeric_limits<size_t>::max()) {
        return true;
    }
    size_t max_trav_len = 0;
    for (const Traversal& walk : walks) {
        // Count the interior nodes only.
        size_t trav_len = 0;
        for (size_t i = 1; i + 1 < walk.size(); ++i) {
            trav_len += graph.get_length(walk[i]);
        }
        max_trav_len = max(max_trav_len, trav_len);
        if (max_trav_len > allele_length_range.second) {
            return false;
        }
    }
    return max_trav_len >= allele_length_range.first;
}

/**
 * Choose the reference walk of a site with no reference path of its own, adding it to `walks` if
 * needed.
 *
 * The first walk of the first nonempty parent set stands in. If no set names one, use `walks[0]`,
 * with parent sets, or none, without them. That leaves a site of off-reference nested calling with
 * no reference walk, since `walks[0]` is the flow finder's best-supported walk, and making it REF
 * would bias the genotyper's tie-break toward it.
 *
 * @return The index in `walks` of the reference walk, or -1 for none.
 */
static int stand_in_reference(vector<Traversal>& travs,
                              const vector<WalkSet>* parent_child_walk_sets) {
    if (parent_child_walk_sets == nullptr) {
        return -1;
    }
    int ref_trav_idx = -1;
    auto first_set = std::ranges::find_if(*parent_child_walk_sets,
                                          [](const WalkSet& set) { return !set.empty(); });
    if (first_set != parent_child_walk_sets->end()) {
        const Traversal& first_trav = first_set->front();
        auto found = std::find(travs.begin(), travs.end(), first_trav);
        if (found != travs.end()) {
            ref_trav_idx = (int)(found - travs.begin());
        } else if (!first_trav.empty()) {
            ref_trav_idx = travs.size();
            travs.push_back(first_trav);
        }
    }
    if (ref_trav_idx < 0) {
        ref_trav_idx = travs.empty() ? -1 : 0;
    }
    return ref_trav_idx;
}

bool CandidateFinder::find(const handle_t& given_start, const handle_t& given_end,
                           const string& parent_ref_path_name,
                           pair<size_t, size_t> parent_ref_interval,
                           const vector<WalkSet>* parent_child_walk_sets, bool no_reference,
                           WalkSite& site) const {
    // The bounds can be turned round to run forward along the reference.
    site.start = given_start;
    site.end = given_end;

    if (graph.get_id(given_start) == graph.get_id(given_end)) {
        // A site of one node cannot be called.
        return false;
    }
    if (too_big(given_start, given_end)) {
        return false;
    }

    // A VCF record needs a reference path through the site.
    const vector<string> common_names = common_reference_paths(given_start, given_end);
    // A site with no reference path takes its parent's path and interval, with parent walk sets or
    // with `no_reference`. Otherwise get_ref_interval would be called with the parent's path,
    // which does not visit this site's boundary nodes, and would assert.
    site.use_parent_interval = common_names.empty();
    if (site.use_parent_interval
        && ((parent_child_walk_sets == nullptr && !no_reference) || parent_ref_path_name.empty())) {
        return false;
    }

    Traversal ref_walk;
    if (site.use_parent_interval) {
        site.ref_path_name = parent_ref_path_name;
        site.ref_interval = make_tuple(parent_ref_interval.first, parent_ref_interval.second, false,
                                       step_handle_t(), step_handle_t());
    } else {
        site.ref_path_name = preferred_reference(common_names);
        site.ref_interval = get_ref_interval(graph, site.start, site.end, site.ref_path_name);
        if (get<0>(site.ref_interval) == -1) {
            // The reference path's orientation conflicts with the site's.
            return false;
        }
        if (get<2>(site.ref_interval)) {
            // Turn the site round to run forward along the reference.
            const handle_t old_start = site.start;
            site.start = graph.flip(site.end);
            site.end = graph.flip(old_start);
            site.ref_interval = get_ref_interval(graph, site.start, site.end, site.ref_path_name);
        }
        if (!reference_walk(site.start, site.end, site.ref_interval, ref_walk)) {
            return false;
        }
    }

    vector<Traversal>& walks = site.walks;
    walks = traversal_finder.find_traversals(site.start, site.end);
    if (walks.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any "
                "traversals through the site "
             << describe_site(graph, given_start, given_end) << endl;
        return false;
    }

    // The length range can keep a giant site from being resolved.
    if (!in_allele_length_range(walks)) {
        return false;
    }

    if (site.use_parent_interval) {
        site.ref_trav_idx = stand_in_reference(walks, parent_child_walk_sets);
    } else {
        // Add the reference walk if the finder did not find it.
        auto found = std::find(walks.begin(), walks.end(), ref_walk);
        site.ref_trav_idx = (int)(found - walks.begin());
        if (found == walks.end()) {
            walks.push_back(ref_walk);
        }
    }
    return true;
}

bool CandidateFinder::find(const Snarl& given_snarl, const string& parent_ref_path_name,
                           pair<size_t, size_t> parent_ref_interval,
                           const vector<TraversalSet>* parent_child_trav_sets, bool no_reference,
                           Site& site) const {
    site.snarl = given_snarl;
    if (!graph.has_node(given_snarl.start().node_id()) || !graph.has_node(given_snarl.end().node_id())) {
        // A site outside the graph cannot be called.
        return false;
    }
    vector<WalkSet> walk_sets;
    if (parent_child_trav_sets != nullptr) {
        for (const TraversalSet& set : *parent_child_trav_sets) {
            WalkSet& walks = walk_sets.emplace_back();
            for (const SnarlTraversal& trav : set) {
                walks.push_back(walk_of_traversal(graph, trav));
            }
        }
    }
    WalkSite walk_site;
    const bool found = find(
        graph.get_handle(given_snarl.start().node_id(), given_snarl.start().backward()),
        graph.get_handle(given_snarl.end().node_id(), given_snarl.end().backward()),
        parent_ref_path_name, parent_ref_interval,
        parent_child_trav_sets != nullptr ? &walk_sets : nullptr, no_reference, walk_site);
    site.snarl.mutable_start()->set_node_id(graph.get_id(walk_site.start));
    site.snarl.mutable_start()->set_backward(graph.get_is_reverse(walk_site.start));
    site.snarl.mutable_end()->set_node_id(graph.get_id(walk_site.end));
    site.snarl.mutable_end()->set_backward(graph.get_is_reverse(walk_site.end));
    site.ref_path_name = walk_site.ref_path_name;
    site.ref_interval = walk_site.ref_interval;
    site.use_parent_interval = walk_site.use_parent_interval;
    site.travs.clear();
    for (const Traversal& walk : walk_site.walks) {
        site.travs.push_back(traversal_of_walk(graph, walk));
    }
    site.ref_trav_idx = walk_site.ref_trav_idx;
    return found;
}

WalkSet CandidateFinder::find_child_walk_set(const Traversal& parent_walk,
                                             const handle_t& child_start,
                                             const handle_t& child_end) const {
    // The parent walk crosses the child if it visits both of the child's boundary nodes.
    const nid_t child_start_id = graph.get_id(child_start);
    const nid_t child_end_id = graph.get_id(child_end);
    bool found_start = false;
    bool found_end = false;
    for (const handle_t& visit : parent_walk) {
        found_start = found_start || graph.get_id(visit) == child_start_id;
        found_end = found_end || graph.get_id(visit) == child_end_id;
    }
    if (!found_start || !found_end) {
        return WalkSet();
    }
    return traversal_finder.find_traversals(child_start, child_end);
}

TraversalSet CandidateFinder::find_child_traversal_set(const SnarlTraversal& parent_trav,
                                                        const Snarl& child) const {
    TraversalSet set;
    for (const Traversal& walk : find_child_walk_set(
             walk_of_traversal(graph, parent_trav),
             graph.get_handle(child.start().node_id(), child.start().backward()),
             graph.get_handle(child.end().node_id(), child.end().backward()))) {
        set.push_back(traversal_of_walk(graph, walk));
    }
    return set;
}

}
