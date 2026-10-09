#include <algorithm>
#include <limits>

#include "candidate_finder.hpp"
#include "graph_caller.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg {

CandidateFinder::CandidateFinder(const PathPositionHandleGraph& graph, SnarlManager& snarl_manager,
                                 TraversalFinder& traversal_finder,
                                 const TraversalSupportFinder& support_finder,
                                 const unordered_set<string>& ref_path_set,
                                 const pair<size_t, size_t>& allele_length_range) :
    graph(graph), snarl_manager(snarl_manager), traversal_finder(traversal_finder),
    support_finder(support_finder), ref_path_set(ref_path_set),
    allele_length_range(allele_length_range) {
}

/// Where `snarl` lies on `ref_path_name`, as `get_ref_interval` gives it for the snarl's bounds.
static tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval_of(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) {
    return get_ref_interval(graph, graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                            graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
                            ref_path_name);
}

bool CandidateFinder::find(const Snarl& managed_snarl, const string& parent_ref_path_name,
                           pair<size_t, size_t> parent_ref_interval,
                           const ChildTraversalSets* parent_child_trav_sets, bool no_reference,
                           Site& site) const {
    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash.
    site.snarl = managed_snarl;
    Snarl& snarl = site.snarl;

#ifdef debug
    cerr << "call_snarl_internal on " << pb2json(snarl) << " with parent_ref_path=" << parent_ref_path_name
         << " parent_child_trav_sets=" << (parent_child_trav_sets ? "provided" : "null") << endl;
#endif

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }

    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    //
    // Only the flow traversal finder uses greedy_avg_flow, so the sum is computed only when there
    // is one.
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.deep_contents(&snarl, graph, false);
        if (snarl_contents.second.size() > max_snarl_edges) {
            // size cap needed as non-nested FlowCaller doesn't handle large snarls
            return false;
        }
        if (flow_trav_finder != nullptr) {
            size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
            size_t length = 0;
            for (auto i = snarl_contents.first.begin();
                 i != snarl_contents.first.end() && length < len_threshold; ++i) {
                length += graph.get_length(graph.get_handle(*i));
            }
            greedy_avg_flow = length > len_threshold;
        }
    }
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    if (common_names.empty()) {
        // No reference path through snarl
        // If we have parent context, we can still process using parent's ref path
        // This test and the use_parent_interval test below must agree: otherwise get_ref_interval
        // would be called with the parent's reference path, which does not visit this snarl's
        // boundary nodes, and would assert.
        if ((parent_child_trav_sets == nullptr && !no_reference)
            || parent_ref_path_name.empty()) {
#ifdef debug
            cerr << "  -> returning false: no common ref path and no parent context" << endl;
#endif
            return false;
        }
#ifdef debug
        cerr << "  -> using parent ref path: " << parent_ref_path_name << endl;
#endif
    }

    // Use parent's ref path if no direct path, otherwise prefer base reference over gref paths
    string& ref_path_name = site.ref_path_name;
    if (common_names.empty()) {
        ref_path_name = parent_ref_path_name;
    } else {
        // Prefer base reference paths over derived gref paths.  Test the whole gref
        // namespace, not just the fragment suffix: a gref copy of the reference sorts
        // before the path it was copied from (gref_x < x).
        // common_names is sorted, so we iterate to find first non-gref path
        ref_path_name = common_names.front();  // default to first (lexicographically smallest)
        for (const string& name : common_names) {
            if (!GrefCover::is_gref_derived(name)) {
                ref_path_name = name;
                break;
            }
        }
    }

    // find the reference traversal and coordinates using the path position graph interface
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t>& ref_interval = site.ref_interval;
    bool& use_parent_interval = site.use_parent_interval;
    use_parent_interval = false;

    if (common_names.empty()) {
        // No direct reference path - use parent's interval and traversals directly
        ref_interval = make_tuple(parent_ref_interval.first, parent_ref_interval.second, false, step_handle_t(), step_handle_t());
        use_parent_interval = true;
    } else {
        ref_interval = ref_interval_of(graph, snarl, ref_path_name);
        if (get<0>(ref_interval) == -1) {
            // could not find reference path interval consistent with snarl due to orientation conflict
            return false;
        }
        if (get<2>(ref_interval) == true) {
            // calling code assumes snarl forward on reference
            flip_snarl(snarl);
            ref_interval = ref_interval_of(graph, snarl, ref_path_name);
        }
    }

    SnarlTraversal ref_trav;

    if (!use_parent_interval) {
        // Build reference traversal from path steps
        step_handle_t cur_step = get<3>(ref_interval);
        step_handle_t last_step = get<4>(ref_interval);
        if (get<2>(ref_interval)) {
            std::swap(cur_step, last_step);
        }
        bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));

        while (true) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(cur_handle));
            visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
            if (graph.get_id(cur_handle) == snarl.end().node_id()) {
                break;
            } else if (get<2>(ref_interval) == true) {
                if (!graph.has_previous_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_previous_step(cur_step);
            } else {
                if (!graph.has_next_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                    return false;
                }
                cur_step = graph.get_next_step(cur_step);
            }
            // todo: we can compute flow at the same time
        }
        assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());
    }
    // If use_parent_interval, ref_trav stays empty - we'll use first parent traversal as pseudo-reference

    vector<SnarlTraversal>& travs = site.travs;
    travs.clear();
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow);
        travs = std::move(weighted_travs.first);
    } else {
        // find the traversals using the generic interface
        travs = traversal_finder.find_traversals(snarl);
    }

    if (travs.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any traversals through snarl " << pb2json(managed_snarl) << endl;
        return false;
    }
#ifdef debug
    cerr << "  found " << travs.size() << " traversals, use_parent_interval=" << use_parent_interval << endl;
#endif

    // optional traversal length clamp can, ex, avoid trying to resolve a giant snarl    
    if (allele_length_range.first > 0 || allele_length_range.second < numeric_limits<size_t>::max()) {
        size_t max_trav_len = 0;
        for (const SnarlTraversal & trav : travs) {
            size_t trav_len = 0;
            for (size_t i = 1; i < trav.visit_size() - 1; ++i) {
                trav_len += graph.get_length(graph.get_handle(trav.visit(i).node_id()));
            }
            max_trav_len = max(max_trav_len, trav_len);
            if (max_trav_len > allele_length_range.second) {
                return false;
            }
        }
        if (max_trav_len < allele_length_range.first) {
            return false;
        }
    }

    // find the reference traversal in the list of results from the traversal finder
    int& ref_trav_idx = site.ref_trav_idx;
    ref_trav_idx = -1;

    if (use_parent_interval) {
        // No direct reference path - use first traversal from first non-empty set as pseudo-reference
        if (parent_child_trav_sets != nullptr) {
            for (const auto& tset : *parent_child_trav_sets) {
                if (!tset.empty()) {
                    const SnarlTraversal& first_trav = tset[0];
                    for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
                        if (travs[i] == first_trav) {
                            ref_trav_idx = i;
                        }
                    }
                    if (ref_trav_idx < 0 && first_trav.visit_size() > 0) {
                        ref_trav_idx = travs.size();
                        travs.push_back(first_trav);
                    }
                    break;
                }
            }
        }
        // Left at -1 where no parent traversal set named a reference: travs[0] is the flow finder's
        // best-supported traversal, and making it REF would bias the genotyper's tie-break toward
        // it. The genotyper checks ref_trav_idx >= 0 before using it.
        if (ref_trav_idx < 0 && parent_child_trav_sets != nullptr) {
            ref_trav_idx = travs.empty() ? -1 : 0;
        }
    } else {
        for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
            // todo: is there a way to speed this up?
            if (travs[i] == ref_trav) {
                ref_trav_idx = i;
            }
        }

        if (ref_trav_idx == -1) {
            ref_trav_idx = travs.size();
            // we didn't get the reference traversal from the finder, so we add it here
            travs.push_back(ref_trav);
        }
    }
    return true;
}
TraversalSet CandidateFinder::find_child_traversal_set(const SnarlTraversal& parent_trav,
                                                        const Snarl& child) const {
    TraversalSet result;

    // First, check if the parent traversal goes through this child snarl
    // by finding the child's start and end nodes in the parent
    nid_t child_start_id = child.start().node_id();
    nid_t child_end_id = child.end().node_id();
    bool found_start = false, found_end = false;

    for (int i = 0; i < parent_trav.visit_size(); ++i) {
        nid_t visit_id = parent_trav.visit(i).node_id();
        if (visit_id == child_start_id) found_start = true;
        if (visit_id == child_end_id) found_end = true;
    }

    // If parent doesn't traverse the child, return empty set (star allele case)
    if (!found_start || !found_end) {
        return result;
    }

    // Use the traversal finder to enumerate all traversals through the child
    FlowTraversalFinder* flow_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_finder != nullptr) {
        auto weighted_travs = flow_finder->find_weighted_traversals(child, false);
        result = std::move(weighted_travs.first);
    } else {
        result = traversal_finder.find_traversals(child);
    }

    return result;
}

}
