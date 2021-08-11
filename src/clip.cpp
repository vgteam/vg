#include "clip.hpp"
#include "traversal_finder.hpp"
#include <unordered_map>
#include <IntervalTree.h>

//#define debug

namespace vg {

using namespace std;

// find the snarl's spanning interval on every reference path traversal through it
// as soon as one of these intervals is found that is contained within an interval in the input index
// then return it (or nullptr if none found)
// also return the snarl's interval (as pair of offsets) in the path
// this logic is mostly lifted from deconstructor which does the same thing to get vcf coordinates.
static pair<const Region*, pair<step_handle_t, step_handle_t>> get_containing_region(PathPositionHandleGraph* graph,
                                                                                     PathTraversalFinder& trav_finder,
                                                                                     const Snarl* snarl,
                                                                                     unordered_map<string, IntervalTree<int64_t, const Region*>>& contig_to_interval_tree,
                                                                                     bool include_endpoints) {

    // every path through the snarl
    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > travs = trav_finder.find_path_traversals(*snarl);

    // check each one against the interval tree
    for (size_t i = 0; i < travs.first.size(); ++i) {
        auto& step_pair = travs.second[i];
        auto& ref_trav = travs.first[i];

        path_handle_t path_handle = graph->get_path_handle_of_step(step_pair.first);
        string path_name = graph->get_path_name(path_handle);
        int64_t path_offset = 0;
        auto sp_parse = Paths::parse_subpath_name(path_name);
        if (get<0>(sp_parse)) {
            path_name = get<1>(sp_parse);
            path_offset = get<2>(sp_parse);
        }

        IntervalTree<int64_t, const Region*>& interval_tree = contig_to_interval_tree.at(path_name);

        // first_path_pos computation copied from deconstructor.cpp (it does not include the start node)
        step_handle_t start_step = step_pair.first;
        step_handle_t end_step = step_pair.second;
        handle_t start_handle = graph->get_handle_of_step(start_step);
        handle_t end_handle = graph->get_handle_of_step(end_step);
        size_t start_pos = graph->get_position_of_step(start_step);
        size_t end_pos = graph->get_position_of_step(end_step);
        bool use_start = start_pos < end_pos;
        handle_t first_path_handle = use_start ? start_handle : end_handle;
        int64_t first_path_pos = use_start ? start_pos : end_pos;
        // Get the first visit of our snarl traversal
        const Visit& first_trav_visit = use_start ? ref_trav.visit(0) : ref_trav.visit(ref_trav.visit_size() - 1);
        if ((use_start && first_trav_visit.backward() == graph->get_is_reverse(first_path_handle)) ||
            (!use_start && first_trav_visit.backward() != graph->get_is_reverse(first_path_handle))) {
            // Our path and traversal have consistent orientation.  leave off the end of the start node going forward
            first_path_pos += graph->get_length(first_path_handle);
        }

        size_t length_from_start = 0;
        for (size_t j = 1; j < ref_trav.visit_size() - 1; ++j) {
            length_from_start += graph->get_length(graph->get_handle(ref_trav.visit(j).node_id()));
        }

        if (include_endpoints) {
            first_path_pos -= graph->get_length(first_path_handle);
            length_from_start += graph->get_length(graph->get_handle(ref_trav.visit(ref_trav.visit_size() - 1).node_id()));
        }
        int64_t last_path_pos = length_from_start == 0 ? first_path_pos : first_path_pos + length_from_start - 1;
        auto overlapping_intervals = interval_tree.findOverlapping(first_path_pos, last_path_pos);
        for (auto& interval : overlapping_intervals) {
            if (interval.start <= first_path_pos && interval.stop >= last_path_pos) {
                return make_pair(interval.value, make_pair(start_step, end_step));
            }
        }
                                                                                                                     
    }
    return make_pair(nullptr, make_pair(step_handle_t(), step_handle_t()));
}

void visit_contained_snarls(PathPositionHandleGraph* graph, const vector<Region>& regions, SnarlManager& snarl_manager,
                            bool include_endpoints,
                            function<void(const Snarl*, const step_handle_t&, const step_handle_t&, const Region*)> visit_fn) {

    // make an interval tree of regions for each contig
    unordered_map<string, vector<IntervalTree<int64_t, const Region*>::interval>> region_intervals;
    for (const Region & region : regions) {
        vector<IntervalTree<int64_t, const Region*>::interval>& intervals = region_intervals[region.seq];
        intervals.push_back(IntervalTree<int64_t, const Region*>::interval(region.start, region.end, &region));
    }
    unordered_map<string, IntervalTree<int64_t, const Region*>> contig_to_interval_tree;
    for (auto seq_intervals : region_intervals) {
        IntervalTree<int64_t, const Region*> interval_tree(std::move(seq_intervals.second));
        contig_to_interval_tree[seq_intervals.first] = std::move(interval_tree);
    }
    region_intervals.clear();
    
    // make a path traversal finder for all affected reference paths taking into account
    // subpaths in the graph
    unordered_set<string> path_name_set;
    for (const Region& region : regions) {
        path_name_set.insert(region.seq);
    }
    unordered_set<string> graph_path_name_set;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            auto sp_parse = Paths::parse_subpath_name(path_name);
            if (get<0>(sp_parse)) {
                path_name = get<1>(sp_parse);
            }
            if (path_name_set.count(path_name)) {
                graph_path_name_set.insert(path_name);
            }
        });
    vector<string> path_names;
    for (const string& path_name : path_name_set) {
        path_names.push_back(path_name);
    }
    path_name_set.clear();
    graph_path_name_set.clear();
    PathTraversalFinder trav_finder(*graph, snarl_manager, path_names);
    
    // Do the top-level snarls, the recurse as needed (framework copied from deconstructor.cpp)
    snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            vector<const Snarl*> todo(1, snarl);
            vector<const Snarl*> next;
            while (!todo.empty()) {
                for (auto next_snarl : todo) {
                    auto containing_region_info = get_containing_region(graph, trav_finder, next_snarl, contig_to_interval_tree, include_endpoints);
                    if (containing_region_info.first != nullptr) {
                        visit_fn(next_snarl, containing_region_info.second.first, containing_region_info.second.second, containing_region_info.first);
                    } else {
                        const vector<const Snarl*>& children = snarl_manager.children_of(next_snarl);
                        next.insert(next.end(), children.begin(), children.end());
                    }
                }
                swap(todo, next);
                next.clear();
            }
        });
}

// note: end_step is after the subpath
//       end_offset is also one-past the interval
static path_handle_t create_path_fragment(MutablePathMutableHandleGraph* graph, const string& base_name, step_handle_t first_step,
                                          step_handle_t end_step, int64_t start_offset, int64_t end_offset) {
    assert(end_offset > start_offset);
    string subpath_name = Paths::make_subpath_name(base_name, start_offset, end_offset - 1);
#ifdef debug
    cerr << "making fragment " << subpath_name << endl;
#endif
    path_handle_t subpath_handle = graph->create_path_handle(subpath_name);
    for (step_handle_t step = first_step; step != end_step; step = graph->get_next_step(step)) {
        graph->append_step(subpath_handle, graph->get_handle_of_step(step));
    }
    return subpath_handle;
}

void clip_contained_snarls(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                           SnarlManager& snarl_manager, bool include_endpoints, bool verbose) {

    // find all nodes in the snarl that are not on the reference interval (reference path name from containing interval)
    unordered_set<nid_t> to_delete;

    // just for logging
    unordered_map<string, size_t> clip_counts;
    
    visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, const step_handle_t& start_step, const step_handle_t& end_step, const Region* containing_region) {

#ifdef debug
            cerr << "Clipping snarl " << pb2json(*snarl) << " because it lies in region "
                 << containing_region->seq << ":" << containing_region->start << "-" << containing_region->end << endl;
#endif

            unordered_set<nid_t> whitelist;
            step_handle_t past_end_step = pp_graph->get_next_step(end_step);            
            for (step_handle_t step = start_step ; step != past_end_step; step = graph->get_next_step(step)) {
                whitelist.insert(pp_graph->get_id(pp_graph->get_handle_of_step(step)));
            }

            pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
            for (id_t node_id : contents.first) {
                if (!whitelist.count(node_id)) {
                    to_delete.insert(node_id);
                    ++clip_counts[containing_region->seq];
                }
            }
            
        });

    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Removing " << kv.second << " nodes due to intervals on path " << kv.first << endl;
        }
        clip_counts.clear();
    }
    
    // chop the paths    
    vector<path_handle_t> path_handles;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });    
    for (path_handle_t& path_handle : path_handles) {

        string path_name = graph->get_path_name(path_handle);
#ifdef debug
        cerr << "processing path " << path_name << endl;
#endif
        int64_t path_offset = 0;
        auto sp_parse = Paths::parse_subpath_name(path_name);
        if (get<0>(sp_parse)) {
            path_name = get<1>(sp_parse);
            path_offset = get<2>(sp_parse);
        }

        int64_t cur_step_offset = path_offset;
        step_handle_t start_step = graph->path_begin(path_handle);
        int64_t start_step_offset = cur_step_offset;
        step_handle_t prev_step;
        step_handle_t end_step = graph->path_end(path_handle);
        bool in_path = !to_delete.count(graph->get_id(graph->get_handle_of_step(start_step)));
        bool was_chopped = false;
        for (step_handle_t cur_step = start_step; cur_step != end_step; cur_step = graph->get_next_step(cur_step)) {
            handle_t cur_handle = graph->get_handle_of_step(cur_step);
            nid_t cur_id = graph->get_id(cur_handle);
            bool step_deleted = to_delete.count(cur_id);
            if (in_path && step_deleted && cur_step_offset > start_step_offset) {
                // we hit a deleted node: make a path fragment for eveything to it
                create_path_fragment(graph, path_name, start_step, cur_step, start_step_offset, cur_step_offset);
                ++clip_counts[path_name];
            }
            if (!in_path && !step_deleted) {
                // we hit the first undeleted node after a run of deleteions, start a new fragment
                start_step_offset = cur_step_offset;
                start_step = cur_step;
            }
            in_path = !step_deleted;
            cur_step_offset += graph->get_length(cur_handle);
            prev_step = cur_step;
            was_chopped = was_chopped || step_deleted;
        }

        if (was_chopped && in_path && cur_step_offset > start_step_offset) {
            // get that last fragment
            create_path_fragment(graph, path_name, start_step, graph->path_end(path_handle), start_step_offset, cur_step_offset);
            ++clip_counts[path_name];
        }

        if (was_chopped) {
            graph->destroy_path(path_handle);                
        }
    }

    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path" << kv.first << endl;
        }
        clip_counts.clear();
    }

    // finally, delete the nodes
    DeletableHandleGraph* del_graph = dynamic_cast<DeletableHandleGraph*>(graph);
    for (nid_t node_id : to_delete) {
        handle_t handle = graph->get_handle(node_id);
        assert(graph->steps_of_handle(handle).empty());
        dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(handle);
    }    
}




}
