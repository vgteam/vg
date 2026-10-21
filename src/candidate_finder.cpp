#include <algorithm>
#include <limits>

#include "candidate_finder.hpp"
#include "graph_caller.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg {

CandidateFinder::CandidateFinder(const PathPositionHandleGraph& graph,
                                 TraversalFinder& traversal_finder,
                                 const unordered_set<string>& ref_path_set,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph), traversal_finder(traversal_finder),
    usable_path([&ref_path_set](const string& name) {
        return ref_path_set.empty() || ref_path_set.contains(name);
    }),
    allele_length_range(allele_length_range) {
}

CandidateFinder::CandidateFinder(const PathPositionHandleGraph& graph,
                                 TraversalFinder& traversal_finder,
                                 const multipass::ReferencePathTable& ref_paths,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph), traversal_finder(traversal_finder),
    usable_path([&ref_paths](const string& name) {
        return ref_paths.empty() || ref_paths.contains(name);
    }),
    allele_length_range(allele_length_range) {
}

/// Find where `snarl` lies on `ref_path_name`, as `get_ref_interval` gives it for the bounds of the
/// snarl.
static tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval_of(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) {
    return get_ref_interval(graph, graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                            graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
                            ref_path_name);
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

bool CandidateFinder::reference_traversal(const Snarl& snarl,
                                          const tuple<int64_t, int64_t, bool, step_handle_t,
                                                      step_handle_t>& ref_interval,
                                          const Snarl& given_snarl, SnarlTraversal& ref_trav) const {
    step_handle_t cur_step = get<3>(ref_interval);
    step_handle_t last_step = get<4>(ref_interval);
    if (get<2>(ref_interval)) {
        std::swap(cur_step, last_step);
    }
    const bool start_backwards =
        snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));

    while (true) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        Visit* visit = ref_trav.add_visit();
        visit->set_node_id(graph.get_id(cur_handle));
        visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle)
                                            : graph.get_is_reverse(cur_handle));
        if (graph.get_id(cur_handle) == snarl.end().node_id()) {
            break;
        }
        const bool path_backward = get<2>(ref_interval);
        if (path_backward ? !graph.has_previous_step(cur_step) : !graph.has_next_step(cur_step)) {
            cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace "
                    "reference path through snarl "
                 << pb2json(given_snarl) << endl;
            return false;
        }
        cur_step = path_backward ? graph.get_previous_step(cur_step) : graph.get_next_step(cur_step);
    }
    assert(ref_trav.visit(0) == snarl.start()
           && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());
    return true;
}

bool CandidateFinder::in_allele_length_range(const vector<SnarlTraversal>& travs) const {
    if (allele_length_range.first == 0
        && allele_length_range.second == numeric_limits<size_t>::max()) {
        return true;
    }
    size_t max_trav_len = 0;
    for (const SnarlTraversal& trav : travs) {
        // Count the interior nodes only.
        size_t trav_len = 0;
        for (size_t i = 1; i < trav.visit_size() - 1; ++i) {
            trav_len += graph.get_length(graph.get_handle(trav.visit(i).node_id()));
        }
        max_trav_len = max(max_trav_len, trav_len);
        if (max_trav_len > allele_length_range.second) {
            return false;
        }
    }
    return max_trav_len >= allele_length_range.first;
}

/**
 * Choose the reference traversal of a site with no reference path of its own, adding it to
 * `travs` if needed.
 *
 * The first traversal of the first nonempty parent set stands in. If no set names one, use
 * `travs[0]`, with parent sets, or none, without them. That leaves a site of off-reference nested
 * calling with no reference traversal, since `travs[0]` is the flow finder's best-supported
 * traversal, and making it REF would bias the genotyper's tie-break toward it.
 *
 * @return The index in `travs` of the reference traversal, or -1 for none.
 */
static int stand_in_reference(vector<SnarlTraversal>& travs,
                              const ChildTraversalSets* parent_child_trav_sets) {
    if (parent_child_trav_sets == nullptr) {
        return -1;
    }
    int ref_trav_idx = -1;
    auto first_set = std::ranges::find_if(*parent_child_trav_sets,
                                          [](const TraversalSet& tset) { return !tset.empty(); });
    if (first_set != parent_child_trav_sets->end()) {
        const SnarlTraversal& first_trav = first_set->front();
        auto found = std::find(travs.begin(), travs.end(), first_trav);
        if (found != travs.end()) {
            ref_trav_idx = (int)(found - travs.begin());
        } else if (first_trav.visit_size() > 0) {
            ref_trav_idx = travs.size();
            travs.push_back(first_trav);
        }
    }
    if (ref_trav_idx < 0) {
        ref_trav_idx = travs.empty() ? -1 : 0;
    }
    return ref_trav_idx;
}

bool CandidateFinder::find(const Snarl& given_snarl, const string& parent_ref_path_name,
                           pair<size_t, size_t> parent_ref_interval,
                           const ChildTraversalSets* parent_child_trav_sets, bool no_reference,
                           Site& site) const {
    // Work on a copy, which can be flipped.
    site.snarl = given_snarl;
    Snarl& snarl = site.snarl;

#ifdef debug
    cerr << "call_snarl_internal on " << pb2json(snarl) << " with parent_ref_path=" << parent_ref_path_name
         << " parent_child_trav_sets=" << (parent_child_trav_sets ? "provided" : "null") << endl;
#endif

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // A site of one node, or outside the graph, cannot be called.
        return false;
    }

    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    if (too_big(start_handle, end_handle)) {
        return false;
    }

    // A VCF record needs a reference path through the site.
    const vector<string> common_names = common_reference_paths(start_handle, end_handle);
    // A site with no reference path takes its parent's path and interval, with parent traversal
    // sets or with `no_reference`. Otherwise get_ref_interval would be called with the parent's
    // path, which does not visit this site's boundary nodes, and would assert.
    site.use_parent_interval = common_names.empty();
    if (site.use_parent_interval
        && ((parent_child_trav_sets == nullptr && !no_reference) || parent_ref_path_name.empty())) {
#ifdef debug
        cerr << "  -> returning false: no common ref path and no parent context" << endl;
#endif
        return false;
    }

    SnarlTraversal ref_trav;
    if (site.use_parent_interval) {
        site.ref_path_name = parent_ref_path_name;
        site.ref_interval = make_tuple(parent_ref_interval.first, parent_ref_interval.second, false,
                                       step_handle_t(), step_handle_t());
    } else {
        site.ref_path_name = preferred_reference(common_names);
        site.ref_interval = ref_interval_of(graph, snarl, site.ref_path_name);
        if (get<0>(site.ref_interval) == -1) {
            // The reference path's orientation conflicts with the site's.
            return false;
        }
        if (get<2>(site.ref_interval)) {
            // Orient the site forward along the reference.
            flip_snarl(snarl);
            site.ref_interval = ref_interval_of(graph, snarl, site.ref_path_name);
        }
        if (!reference_traversal(snarl, site.ref_interval, given_snarl, ref_trav)) {
            return false;
        }
    }

    vector<SnarlTraversal>& travs = site.travs;
    travs = traversal_finder.find_traversals(snarl);
    if (travs.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any traversals through snarl " << pb2json(given_snarl) << endl;
        return false;
    }
#ifdef debug
    cerr << "  found " << travs.size() << " traversals, use_parent_interval=" << site.use_parent_interval << endl;
#endif

    // The length range can keep a giant site from being resolved.
    if (!in_allele_length_range(travs)) {
        return false;
    }

    if (site.use_parent_interval) {
        site.ref_trav_idx = stand_in_reference(travs, parent_child_trav_sets);
    } else {
        // Add the reference traversal if the finder did not find it.
        auto found = std::find(travs.begin(), travs.end(), ref_trav);
        site.ref_trav_idx = (int)(found - travs.begin());
        if (found == travs.end()) {
            travs.push_back(ref_trav);
        }
    }
    return true;
}

TraversalSet CandidateFinder::find_child_traversal_set(const SnarlTraversal& parent_trav,
                                                        const Snarl& child) const {
    // The parent traversal crosses the child if it visits both of the child's boundary nodes.
    const nid_t child_start_id = child.start().node_id();
    const nid_t child_end_id = child.end().node_id();
    bool found_start = false;
    bool found_end = false;
    for (const Visit& visit : parent_trav.visit()) {
        found_start = found_start || visit.node_id() == child_start_id;
        found_end = found_end || visit.node_id() == child_end_id;
    }
    if (!found_start || !found_end) {
        return TraversalSet();
    }
    return traversal_finder.find_traversals(child);
}

}
