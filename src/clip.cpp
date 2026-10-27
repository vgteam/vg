#include "clip.hpp"
#include "traversal_finder.hpp"
#include <unordered_map>
#include <IntervalTree.h>
#include <structures/rank_pairing_heap.hpp>
#include <structures/union_find.hpp>
#include <BooPHF.h>
#include "bdsg/internal/hash_map.hpp"
#include "bdsg/internal/packed_structs.hpp"

//#define debug

namespace vg {

using namespace std;

// find the snarl's spanning interval on every reference path traversal through it
// as soon as one of these intervals is found that is contained within an interval in the input index
// then return it (or nullptr if none found)
// also return the snarl's interval (as pair of offsets) in the path
// this logic is mostly lifted from deconstructor which does the same thing to get vcf coordinates.
static tuple<const Region*, step_handle_t, step_handle_t, int64_t, int64_t, bool> get_containing_region(const PathPositionHandleGraph* graph,
                                                                                                        PathTraversalFinder& trav_finder,
                                                                                                        const Snarl* snarl,
                                                                                                        unordered_map<string, IntervalTree<int64_t, const Region*>>& contig_to_interval_tree,
                                                                                                        bool include_endpoints) {
    
    // every path through the snarl
    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > travs = trav_finder.find_path_traversals(*snarl);

    // we sort by (region-size, interval-size) to choose the biggest interval from the biggest contig
    // intuition: pggb graph with tons of unplaced contigs and one reference contig -- we want the reference contig
    // then in all the possible traversals of that reference contig, we take the biggest (which should fit better
    // with the other heuristics)
    multimap<pair<int64_t, int64_t>,
             tuple<const Region*, step_handle_t, step_handle_t, int64_t, int64_t, bool>> ranked_intervals;
    
    // check each one against the interval tree
    for (size_t i = 0; i < travs.first.size(); ++i) {
        auto& step_pair = travs.second[i];
        auto& ref_trav = travs.first[i];

        path_handle_t path_handle = graph->get_path_handle_of_step(step_pair.first);
        string path_name = graph->get_path_name(path_handle);
        int64_t path_offset = 0;
        subrange_t subrange;
        path_name = Paths::strip_subrange(path_name, &subrange);
        if (subrange != PathMetadata::NO_SUBRANGE) {
            path_offset = subrange.first;
        }

        if (contig_to_interval_tree.count(path_name)) {

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
            int64_t traversal_interval_length = last_path_pos - first_path_pos + 1;            
            for (auto& interval : overlapping_intervals) {
                if (interval.start <= first_path_pos && interval.stop >= last_path_pos) {
                    int64_t region_interval_length = interval.stop - interval.start + 1;
                    ranked_intervals.insert(make_pair(make_pair(region_interval_length, traversal_interval_length),
                                            make_tuple(interval.value, start_step, end_step, first_path_pos, last_path_pos, !use_start)));
                }
            }
        }
    }

    
    if (!ranked_intervals.empty()) {
        return ranked_intervals.rbegin()->second;
    }
    return make_tuple(nullptr, step_handle_t(), step_handle_t(), -1, -1, false);
}

void visit_contained_snarls(const PathPositionHandleGraph* graph, const vector<Region>& regions, SnarlManager& snarl_manager,
                            bool include_endpoints,
                            function<void(const Snarl*, step_handle_t, step_handle_t, int64_t, int64_t, bool, const Region*)> visit_fn) {

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
            string graph_path_name = graph->get_path_name(path_handle);
            if (path_name_set.count(Paths::strip_subrange(graph_path_name))) {
                graph_path_name_set.insert(graph_path_name);
            }
        });
    vector<string> path_names;
    for (const string& path_name : graph_path_name_set) {
        path_names.push_back(path_name);
    }
    path_name_set.clear();
    graph_path_name_set.clear();
    PathTraversalFinder trav_finder(*graph, path_names);
    
    // Do the top-level snarls, the recurse as needed (framework copied from deconstructor.cpp)
    snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            vector<const Snarl*> todo(1, snarl);
            vector<const Snarl*> next;
            while (!todo.empty()) {
                for (auto next_snarl : todo) {
                    auto containing_region_info = get_containing_region(graph, trav_finder, next_snarl, contig_to_interval_tree, include_endpoints);
                    if (get<0>(containing_region_info)  != nullptr) {
                        visit_fn(next_snarl, get<1>(containing_region_info), get<2>(containing_region_info), get<3>(containing_region_info),
                                 get<4>(containing_region_info), get<5>(containing_region_info), get<0>(containing_region_info));
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
    PathSense sense;
    std::string sample;
    std::string locus;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(base_name, sense, sample, locus, haplotype, phase_block, subrange);
    assert(subrange == PathMetadata::NO_SUBRANGE);
    subrange.first = start_offset;
    subrange.second = end_offset;
    string subpath_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block, subrange);
#ifdef debug
    cerr << "making fragment " << subpath_name << endl;
#endif
    path_handle_t subpath_handle = graph->create_path_handle(subpath_name);
    for (step_handle_t step = first_step; step != end_step; step = graph->get_next_step(step)) {
        graph->append_step(subpath_handle, graph->get_handle_of_step(step));
    }
    return subpath_handle;
}

// note: clip-vg.cpp has a more general version (if ever needed) that can chop nodes on path positions
void delete_nodes_and_chop_paths(MutablePathMutableHandleGraph* graph, const unordered_set<nid_t>& nodes_to_delete,
                                 const unordered_set<edge_t>& edges_to_delete, int64_t min_fragment_len,
                                 unordered_map<string, size_t>* fragments_per_path, bool keep_short_at_deleted_nodes) {
      // chop the paths    
    vector<path_handle_t> path_handles;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });
    vector<path_handle_t> paths_to_destroy;
    for (path_handle_t& path_handle : path_handles) {

        string path_name = graph->get_path_name(path_handle);
#ifdef debug
        cerr << "processing path " << path_name << endl;
#endif
        int64_t path_offset = 0;
        subrange_t subrange;
        path_name = Paths::strip_subrange(path_name, &subrange);
        if (subrange != PathMetadata::NO_SUBRANGE) {
            path_offset = subrange.first;
        }

        int64_t cur_step_offset = path_offset;
        step_handle_t start_step = graph->path_begin(path_handle);
        int64_t start_step_offset = cur_step_offset;
        step_handle_t prev_step;
        step_handle_t end_step = graph->path_end(path_handle);
        bool in_path = !nodes_to_delete.count(graph->get_id(graph->get_handle_of_step(start_step)));
        bool was_chopped = false;
        for (step_handle_t cur_step = start_step; cur_step != end_step; cur_step = graph->get_next_step(cur_step)) {
            handle_t cur_handle = graph->get_handle_of_step(cur_step);
            nid_t cur_id = graph->get_id(cur_handle);
            bool step_deleted = nodes_to_delete.count(cur_id);
            bool edge_deleted = false;
            if (in_path && cur_step_offset > start_step_offset) {
                if (!step_deleted) {
                    edge_deleted = edges_to_delete.count(graph->edge_handle(graph->get_handle_of_step(prev_step), cur_handle));
                }
                if (step_deleted || edge_deleted) {
                    // we hit a deleted node (or edge): make a path fragment for eveything to it
                    if ((step_deleted && keep_short_at_deleted_nodes) || cur_step_offset - start_step_offset >= min_fragment_len) {
                        create_path_fragment(graph, path_name, start_step, cur_step, start_step_offset, cur_step_offset);
                        if (fragments_per_path) {
                            ++(*fragments_per_path)[path_name];
                        }
                    }
                }
                if (edge_deleted) { // need to handle this case by popping off the path now
                    in_path = false;
                }
            }
            if (!in_path && !step_deleted) {
                // we hit the first undeleted node after a run of deleteions, start a new fragment
                start_step_offset = cur_step_offset;
                start_step = cur_step;
            }
            in_path = !step_deleted;
            cur_step_offset += graph->get_length(cur_handle);
            prev_step = cur_step;
            was_chopped = was_chopped || step_deleted || edge_deleted;
        }

        if (was_chopped && in_path && cur_step_offset > start_step_offset) {
            // get that last fragment
            if (cur_step_offset - start_step_offset >= min_fragment_len) {
                create_path_fragment(graph, path_name, start_step, graph->path_end(path_handle), start_step_offset, cur_step_offset);
                if (fragments_per_path) {
                    ++(*fragments_per_path)[path_name];
                }
            }
        }

        if (was_chopped) {
            paths_to_destroy.push_back(path_handle);
        }
    }

    // delete the paths
    graph->destroy_paths(paths_to_destroy);

    DeletableHandleGraph* del_graph = dynamic_cast<DeletableHandleGraph*>(graph);
    // delete the edges
    for (edge_t edge : edges_to_delete) {
        del_graph->destroy_edge(edge);
    }
    
    // finally, delete the nodes    
    for (nid_t node_id : nodes_to_delete) {
        handle_t handle = graph->get_handle(node_id);
        assert(graph->steps_of_handle(handle).empty());
        del_graph->destroy_handle(handle);
    }    

}

// determine if a snarl is complex enough to warrant flattening by measuring some very basic stats
static bool snarl_is_complex(PathPositionHandleGraph* graph, const Snarl* snarl,
                             const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents,
                             const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents_shallow,
                             int64_t ref_interval_length, const Region& region, path_handle_t path_handle,
                             size_t min_nodes, size_t min_edges, size_t min_bases,
                             size_t min_nodes_shallow, size_t min_edges_shallow,
                             double min_avg_degree, double min_reflen, double max_reflen_prop, size_t max_reflen,
                             double& out_avg_degree) {
    
    out_avg_degree = -1.;
    
    // if our snarl is to big vs the reference path, we do not process it
    double ref_prop = (double)ref_interval_length / (double)graph->get_path_length(path_handle);

    if (ref_prop > max_reflen_prop || ref_interval_length > max_reflen || ref_interval_length < min_reflen) {
#ifdef debug
        cerr << "skipping snarl " << pb2json(*snarl) << " with interval length " << ref_interval_length
             << " because its ref_prop of " << region.seq << " is " << ref_prop << " which is greater than " << max_reflen_prop
             << " or its ref length " << ref_interval_length << " is greater than " << max_reflen
             << " or its ref length " << ref_interval_length << " is less than " << min_reflen << endl;
#endif
        return false;
    }

    // check the stats
    bool filter_on = min_nodes > 0 || min_edges > 0 || min_bases > 0 ||
        min_nodes_shallow > 0 || min_edges_shallow > 0 || min_avg_degree > 0.;
    filter_on = filter_on || min_reflen > 0 || max_reflen_prop < 1. || max_reflen < numeric_limits<size_t>::max();
    bool complex_nodes = contents.first.size() > min_nodes;
    bool complex_edges = contents.second.size() > min_edges;
    bool complex_bases = false;
    size_t total_bases = 0;
    for (nid_t node_id : contents.first) {
        total_bases += graph->get_length(graph->get_handle(node_id));
        if (total_bases > min_bases) {
            complex_bases =true;
            break;
        }
    }
    bool complex_nodes_shallow = contents_shallow.first.size() > min_nodes_shallow;
    bool complex_edges_shallow = contents_shallow.second.size() > min_edges_shallow;
    size_t total_degree = 0;
    for (id_t node_id : contents.first) {
        handle_t handle = graph->get_handle(node_id);
        total_degree += graph->get_degree(handle, true) + graph->get_degree(handle, false);
    }
    // degree averaged over node sides to be a bit more intuitive, hence 2X in denominator:
    out_avg_degree = (double)total_degree / (2. *(double)contents.first.size());
    bool complex_degree = out_avg_degree > min_avg_degree;        
    
    return !filter_on || (complex_nodes && complex_edges && complex_bases && complex_nodes_shallow &&
                          complex_edges_shallow && complex_degree);
}

static unordered_set<edge_t> get_net_edges(const HandleGraph* graph, SnarlManager& snarl_manager,
                                           const Snarl* snarl, bool recursive) {
    unordered_set<edge_t> net_edges;

    function<void(const Snarl*)> process_snarl = [&](const Snarl* s) {
        const deque<Chain>& chains = snarl_manager.chains_of(s);

        for (const Chain& chain : chains) {
            const Snarl* first_snarl = chain.front().first;
            if (chain.front().second) {
                snarl_manager.flip(first_snarl);
            }
            handle_t snarl_start = graph->get_handle(first_snarl->start().node_id(), first_snarl->start().backward());
            graph->follow_edges(snarl_start, true, [&](handle_t previous) {
                net_edges.insert(graph->edge_handle(previous, snarl_start));
            });
            if (chain.front().second) {
                snarl_manager.flip(first_snarl);
            }

            const Snarl* last_snarl = chain.back().first;
            if (chain.back().second) {
                snarl_manager.flip(last_snarl);
            }
            handle_t snarl_end = graph->get_handle(last_snarl->end().node_id(), last_snarl->end().backward());
            graph->follow_edges(snarl_end, false, [&](handle_t next) {
                net_edges.insert(graph->edge_handle(snarl_end, next));
            });
            if (chain.back().second) {
                snarl_manager.flip(last_snarl);
            }
        }
    };

    vector<const Snarl*> queue = {snarl};
    while (!queue.empty()) {
        const Snarl* cur_snarl = queue.back();
        queue.pop_back();        
        process_snarl(cur_snarl);    
        if (recursive) {
            const vector<const Snarl*>& children = snarl_manager.children_of(cur_snarl);
            for (const Snarl* child : children) {
                queue.push_back(child);
            }
        }
    }
    
    return net_edges;
}


void clip_contained_snarls(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions, 
                           SnarlManager& snarl_manager, bool include_endpoints, int64_t min_fragment_len,
                           size_t min_nodes, size_t min_edges, size_t min_bases, size_t min_nodes_shallow, size_t min_edges_shallow,
                           double min_avg_degree, size_t min_reflen, double max_reflen_prop, size_t max_reflen, bool only_net_edges,
                           bool only_top_net_edges, bool out_bed, bool verbose) {

    // find all nodes in the snarl that are not on the reference interval (reference path name from containing interval)
    unordered_set<nid_t> nodes_to_delete;

    // and all the edges
    unordered_set<edge_t> edges_to_delete;

    // just for logging
    unordered_map<string, size_t> clip_counts;

    // for making the whitelist
    unordered_set<string> ref_prefixes;
    for (const Region& region : regions) {
        ref_prefixes.insert(region.seq);
    }
    
    visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                                                                                    int64_t start_pos, int64_t end_pos,
                                                                                    bool steps_reversed, const Region* containing_region) {

#ifdef debug
            cerr << "Clipping snarl " << pb2json(*snarl) << " because it lies in region "
                 << containing_region->seq << ":" << containing_region->start << "-" << containing_region->end << endl;
            cerr << "Passed in steps are " << pp_graph->get_id(pp_graph->get_handle_of_step(start_step)) << ":" << pp_graph->get_is_reverse(pp_graph->get_handle_of_step(start_step)) << " - "
                 << pp_graph->get_id(pp_graph->get_handle_of_step(end_step)) << ":" << pp_graph->get_is_reverse(pp_graph->get_handle_of_step(end_step)) << endl;
#endif

            unordered_set<nid_t> whitelist;
            if (steps_reversed) {
                step_handle_t past_end_step = pp_graph->get_previous_step(end_step);            
                for (step_handle_t step = start_step ; step != past_end_step; step = graph->get_previous_step(step)) {
                    whitelist.insert(pp_graph->get_id(pp_graph->get_handle_of_step(step)));
                }
            } else {
                step_handle_t past_end_step = pp_graph->get_next_step(end_step);            
                for (step_handle_t step = start_step ; step != past_end_step; step = graph->get_next_step(step)) {
                    whitelist.insert(pp_graph->get_id(pp_graph->get_handle_of_step(step)));
                }
            }

            edge_t deletion_edge = graph->edge_handle(graph->get_handle(snarl->start().node_id(), snarl->start().backward()),
                                                      graph->get_handle(snarl->end().node_id(), snarl->end().backward()));
            bool deletion_on_whitelist = false;
            // check if the snarl-spanning deletion edge is on a reference path.  if it is, we whitelist it, otherwise it
            // goes.  this gets treated separately as it would not otherwise get deleted by erasing every node in the snarl contents
            if (graph->has_edge(deletion_edge)) {
                graph->for_each_step_on_handle(graph->get_handle(snarl->start().node_id(), snarl->start().backward()), [&](step_handle_t step_handle) {
                    string path_name = graph->get_path_name(graph->get_path_handle_of_step(step_handle));
                    for (const string& ref_prefix : ref_prefixes) {
                        if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                            step_handle_t next_step = graph->get_next_step(step_handle);
                            if (next_step != graph->path_end(graph->get_path_handle_of_step(step_handle)) &&
                                graph->edge_handle(graph->get_handle_of_step(step_handle), graph->get_handle_of_step(next_step)) == deletion_edge) {
                                deletion_on_whitelist = true;
                                return false;
                            }
                            step_handle_t prev_step = graph->get_previous_step(step_handle);
                            if (prev_step != graph->path_front_end(graph->get_path_handle_of_step(step_handle)) &&
                                graph->edge_handle(graph->get_handle_of_step(prev_step), graph->get_handle_of_step(step_handle)) == deletion_edge) {
                                deletion_on_whitelist = true;
                                return false;
                            }
                        }
                    }
                    return true;
                });
            }
            
            size_t ref_interval_length = 0;
            for (nid_t node_id : whitelist) {
                // don't count snarl ends here. todo: should this be an option?
                if (node_id != snarl->start().node_id() && node_id != snarl->end().node_id()) {
                    ref_interval_length += pp_graph->get_length(pp_graph->get_handle(node_id));
                }
            }
            path_handle_t path_handle = pp_graph->get_path_handle_of_step(start_step);
            pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
            pair<unordered_set<id_t>, unordered_set<edge_t> > contents_shallow = snarl_manager.shallow_contents(snarl, *pp_graph, false);
            // add other reference paths to the whitelist to make sure they don't get cut
            if (!ref_prefixes.empty()) {
                for (const id_t& node_id : contents.first) {
                    if (!whitelist.count(node_id)) {
                        graph->for_each_step_on_handle(graph->get_handle(node_id), [&](step_handle_t step_handle) {
                            string path_name = graph->get_path_name(graph->get_path_handle_of_step(step_handle));
                            for (const string& ref_prefix : ref_prefixes) {
                                if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                                    whitelist.insert(node_id);
                                    return false;
                                }
                            }
                            return true;
                        });
                    }
                }
            }
            
            double avg_degree = -1;
            if (snarl_is_complex(pp_graph, snarl, contents, contents_shallow, ref_interval_length, *containing_region, path_handle, min_nodes, min_edges, min_bases,
                                 min_nodes_shallow, min_edges_shallow, min_avg_degree, min_reflen, max_reflen_prop, max_reflen, avg_degree)) {
                if (out_bed) {
                    string snarl_name = (snarl->start().backward() ? "<" : ">") + std::to_string(snarl->start().node_id()) +
                        (snarl->end().backward() ? "<" : ">") + std::to_string(snarl->end().node_id());
                    cout << containing_region->seq << "\t" << start_pos << "\t" << (end_pos + 1) << "\t" << snarl_name << "\t"
                         << contents.first.size() << "\t" << contents.second.size() << "\t"
                         << contents_shallow.first.size() << "\t" << contents_shallow.second.size() << "\t"
                         << avg_degree << "\n";
                } else if (only_net_edges) {
                    unordered_set<edge_t> net_edges = get_net_edges(pp_graph, snarl_manager, snarl, !only_top_net_edges);
                    for (const edge_t& edge : net_edges) {
                        edges_to_delete.insert(edge);
                        assert(pp_graph->has_edge(edge));
                    }
                } else {
                    for (id_t node_id : contents.first) {
                        if (!whitelist.count(node_id)) {
                            nodes_to_delete.insert(node_id);
                            ++clip_counts[containing_region->seq];
                        }
                    }
                    // since we're deleting all alt alleles, the only edge that could be left is a snarl-spanning deletion
                    if (!deletion_on_whitelist && graph->has_edge(deletion_edge)) {
                        edges_to_delete.insert(deletion_edge);
                    }
                }
            }
#ifdef debug
            cerr << "snarl was not deemed complex enough to clip" << endl;
#endif
        });

    if (verbose && !out_bed) {
        if (clip_counts.size() > 1) {
            for (const auto& kv : clip_counts) {
                cerr << "[vg-clip]: Removing " << kv.second << " nodes due to intervals on path " << kv.first << endl;
            }
        }
        cerr << "[vg-clip]: Removing total of " << nodes_to_delete.size() << " nodes and " << edges_to_delete.size() << " snarl-spanning edges from snarls in regions" << endl;
        clip_counts.clear();
    }
    
    // cut out the nodes and chop up paths
    if (!out_bed) {        
        delete_nodes_and_chop_paths(graph, nodes_to_delete, edges_to_delete, min_fragment_len, verbose ? &clip_counts : nullptr);
        if (verbose) {
            for (const auto& kv : clip_counts) {
                cerr << "[vg-clip]: Creating " << kv.second << " fragments from path " << kv.first << endl;
            }
            clip_counts.clear();
        }
    }
}

struct BBEdgeHash {
    uint64_t operator()(const edge_t& edge, uint64_t seed = 0xAAAAAAAA55555555ULL) const {
        uint64_t hsh1 = boomphf::SingleHashFunctor<int64_t>()(as_integer(edge.first), seed);
        uint64_t hsh2 = boomphf::SingleHashFunctor<int64_t>()(as_integer(edge.second), seed);
        // Boost combine for hash values
        return hsh1 ^ (hsh2 + 0x9e3779b9 + (hsh1<<6) + (hsh1>>2));
    }
};

void clip_low_depth_nodes_and_edges_generic(MutablePathMutableHandleGraph* graph,
                                            function<void(function<void(handle_t, const Region*)>)> iterate_handles,
                                            function<void(function<void(edge_t, const Region*)>)> iterate_edges,                                            
                                            int64_t min_depth, const vector<string>& ref_prefixes,
                                            int64_t min_fragment_len, bool verbose) {

    // find all nodes in the snarl that are not on the reference interval (reference path name from containing interval)
    unordered_set<nid_t> to_delete;

    // just for logging
    unordered_map<string, size_t> clip_counts;

    function<bool(const string&)> check_prefixes = [&ref_prefixes] (const string& path_name) {
        for (const string& ref_prefix : ref_prefixes) {
            if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                return true;
            }
        }
        return false;
    };

    function<void(handle_t, const Region*)> visit_handle = [&](handle_t handle, const Region* region) {
        bool on_ref = false;
        size_t depth = 0;
        graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle) {
                ++depth;
                if (depth > min_depth || on_ref) {
                    return false;
                }
                if (!ref_prefixes.empty() || region) {
                    // if we have a region, do exact comparison to it.
                    // otherwise, do a prefix check against ref_prefix
                    string path_name = graph->get_path_name(graph->get_path_handle_of_step(step_handle));
                    if ((region && region->seq == path_name) || (!region && check_prefixes(path_name))) {
                        on_ref = true;
                        return false;
                    }
                }
                return true;
            });
        if (!on_ref && depth < min_depth) {
            to_delete.insert(graph->get_id(handle));
        }
    };

    iterate_handles(visit_handle);
    
    if (verbose) {
        cerr << "[vg-clip]: Removing " << to_delete.size() << " nodes with path coverage less than " << min_depth << endl;
    }

    // now do the edges
    vector<edge_t> edges;
    graph->for_each_edge([&](edge_t edge) {
            edges.push_back(edge);
        });
    size_t edge_count = edges.size();
    boomphf::mphf<edge_t, BBEdgeHash> edge_hash(edge_count, edges, get_thread_count(), 2.0, false, false);
    edges.clear();
    bdsg::PackedVector<> edge_depths;
    edge_depths.resize(edge_count + 1);
    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            bool is_ref_path = check_prefixes(graph->get_path_name(path_handle));
            handle_t prev_handle;
            bool first = true;
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    handle_t handle = graph->get_handle_of_step(step_handle);
                    if (!first) {
                        edge_t edge = graph->edge_handle(prev_handle, handle);
                        size_t edge_rank = edge_hash.lookup(edge);
                        int64_t edge_depth = edge_depths.get(edge_rank);
                        if (edge_depth < min_depth) {
                            if (is_ref_path) {
                                // we never want to remove and edge on a reference path,
                                // so automatically bump such edges past the threshold
                                edge_depths.set(edge_rank, min_depth);
                            } else {
                                edge_depths.set(edge_rank, edge_depth + 1);
                            }
                        }
                    } else {
                        first = false;
                    }
                    prev_handle = handle;
                });
        });

    unordered_set<edge_t> edges_to_delete;
    function<void(edge_t, const Region*)> visit_edge = [&](edge_t edge, const Region* region) {
        size_t edge_rank = edge_hash.lookup(edge);
        if (edge_depths.get(edge_rank) < min_depth) {
            edges_to_delete.insert(edge);
        }
    };

    iterate_edges(visit_edge);

    if (verbose) {
        cerr << "[vg-clip]: Removing " << edges_to_delete.size() << " edges with path coverage less than " << min_depth << endl;
    }

    // cut out the nodes and chop up paths
    delete_nodes_and_chop_paths(graph, to_delete, edges_to_delete, min_fragment_len, verbose ? &clip_counts : nullptr);
        
    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path " << kv.first << endl;
        }
        clip_counts.clear();
    }

    // use the reference path prefix (if given) to clip out components that aren't anchored to it
    // (this would take care of above filter, but we leave that one as it's not dependent on path name)
    if (!ref_prefixes.empty()) {
        size_t removed_node_count = 0;
        size_t removed_component_count = 0;
        vector<unordered_set<nid_t>> components = handlealgs::weakly_connected_components(graph);
        for (auto& component : components) {
            bool ref_anchored = false;
            for (auto ni = component.begin(); !ref_anchored && ni != component.end(); ++ni) {
                vector<step_handle_t> steps = graph->steps_of_handle(graph->get_handle(*ni));
                for (size_t si = 0; !ref_anchored && si < steps.size(); ++si) {
                    string step_path_name = graph->get_path_name(graph->get_path_handle_of_step(steps[si]));
                    if (check_prefixes(step_path_name)) {
                        ref_anchored = true;
                    }
                }
            }
            if (!ref_anchored) {
                ++removed_component_count;
                for (auto node_id : component) {
                    handle_t node_handle = graph->get_handle(node_id);
                    dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(node_handle);
                    ++removed_node_count;
                }
            }
        }
        if (verbose) {
            cerr << "[vg-clip]: Removing " << removed_node_count << " nodes in " << removed_component_count << " disconnected components" << endl;
        }
    }
}

void clip_low_depth_nodes_and_edges(MutablePathMutableHandleGraph* graph, int64_t min_depth, const vector<string>& ref_prefixes,
                          int64_t min_fragment_len, bool verbose) {
    
    function<void(function<void(handle_t, const Region*)>)> iterate_handles = [&] (function<void(handle_t, const Region*)> visit_handle) {
        graph->for_each_handle([&](handle_t handle) {
                visit_handle(handle, nullptr);
            });        
    };

    function<void(function<void(edge_t, const Region*)>)> iterate_edges = [&] (function<void(edge_t, const Region*)> visit_edge) {
        graph->for_each_edge([&](edge_t edge) {
                visit_edge(edge, nullptr);
            });        
    };

    clip_low_depth_nodes_and_edges_generic(graph, iterate_handles, iterate_edges, min_depth, ref_prefixes, min_fragment_len, verbose);
}

void clip_contained_low_depth_nodes_and_edges(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                                              SnarlManager& snarl_manager, bool include_endpoints, int64_t min_depth, int64_t min_fragment_len,
                                              size_t min_nodes, size_t min_edges, size_t min_bases, size_t min_nodes_shallow, size_t min_edges_shallow,
                                              double min_avg_degree, size_t min_reflen, double max_reflen_prop, size_t max_reflen, bool only_net_edges,
                                              bool only_top_net_edges, bool out_bed, bool verbose) {

    vector<pair<Region, pair<vector<handle_t>, vector<edge_t>>>> to_visit;

    visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                                                                                        int64_t start_pos, int64_t end_pos,
                                                                                        bool steps_reversed, const Region* containing_region) {
        unordered_set<nid_t> whitelist;
        if (steps_reversed) {
            step_handle_t past_end_step = pp_graph->get_previous_step(end_step);            
            for (step_handle_t step = start_step ; step != past_end_step; step = graph->get_previous_step(step)) {
                whitelist.insert(pp_graph->get_id(pp_graph->get_handle_of_step(step)));
            }
        } else {
            step_handle_t past_end_step = pp_graph->get_next_step(end_step);            
            for (step_handle_t step = start_step ; step != past_end_step; step = graph->get_next_step(step)) {
                whitelist.insert(pp_graph->get_id(pp_graph->get_handle_of_step(step)));
            }
        }
        size_t ref_interval_length = 0;
        for (nid_t node_id : whitelist) {
            // don't count snarl ends here. todo: should this be an option?
            if (node_id != snarl->start().node_id() && node_id != snarl->end().node_id()) {
                ref_interval_length += pp_graph->get_length(pp_graph->get_handle(node_id));
            }
        }
        path_handle_t path_handle = pp_graph->get_path_handle_of_step(start_step);
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents_shallow = snarl_manager.shallow_contents(snarl, *pp_graph, false);
        double avg_degree = -1;

        if (snarl_is_complex(pp_graph, snarl, contents, contents_shallow, ref_interval_length, *containing_region, path_handle, min_nodes, min_edges, min_bases, 
                             min_nodes_shallow, min_edges_shallow, min_avg_degree, min_reflen, max_reflen_prop, max_reflen, avg_degree)) {
            pair<Region, pair<vector<handle_t>, vector<edge_t>>> region_contents;
            if (containing_region) {
                region_contents.first = *containing_region;
            }            
            if (only_net_edges) {
                unordered_set<edge_t> net_edges = get_net_edges(pp_graph, snarl_manager, snarl, !only_top_net_edges);
                for (const edge_t& edge : net_edges) {
                    region_contents.second.second.push_back(edge);
                    assert(pp_graph->has_edge(edge));
                }
            } else {
                for (id_t node_id : contents.first) {
                    region_contents.second.first.push_back(graph->get_handle(node_id));
                }
                for (const edge_t& edge : contents.second) {
                    region_contents.second.second.push_back(edge);
                }
            }
            to_visit.push_back(region_contents);
        }
    });

    function<void(function<void(handle_t, const Region*)>)> iterate_handles = [&] (function<void(handle_t, const Region*)> visit_handle) {
        for (const auto& region_contents : to_visit) {
            for (const handle_t& handle : region_contents.second.first) {
                visit_handle(handle, &region_contents.first);
            }
        }
    };
    function<void(function<void(edge_t, const Region*)>)> iterate_edges = [&] (function<void(edge_t, const Region*)> visit_edge) {
        for (const auto& region_contents : to_visit) {
            for (const edge_t& edge : region_contents.second.second) {
                visit_edge(edge, &region_contents.first);
            }
        }
    };
    // the edge depths are computed globally, without looking at regions.  as such, they need some notion of reference paths
    // so we shimmy a set in from the regions
    set<string> ref_path_set;
    for (const Region& region : regions) {
        ref_path_set.insert(region.seq);
    }
    vector<string> ref_paths_from_regions(ref_path_set.begin(), ref_path_set.end());

    clip_low_depth_nodes_and_edges_generic(graph, iterate_handles, iterate_edges, min_depth, ref_paths_from_regions, min_fragment_len, verbose);
    
}

// we avoid the path position interface since we only want reference coordinates
static unordered_map<handle_t, vector<int64_t>> make_ref_index(PathHandleGraph* graph, path_handle_t ref_path,
                                                               unordered_set<edge_t>& out_ref_edges) {
    unordered_map<handle_t, vector<int64_t>> handle_to_position;
    int64_t pos = 0;
    handle_t prev_handle;
    bool has_prev = false;
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            int64_t handle_len = graph->get_length(handle);
            handle_to_position[handle].push_back(pos);
            handle_to_position[graph->flip(handle)].push_back(pos + handle_len - 1);
            pos += handle_len;
            if (has_prev) {
                out_ref_edges.insert(graph->edge_handle(prev_handle, handle));
            }
            has_prev = true;
            prev_handle = handle;
        });
    return handle_to_position;
}

// walk context steps out from reference path, flagging each node encountered with its
// minimum and maximum position on the path
static multimap<int64_t, edge_t> find_deletion_candidate_edges(PathHandleGraph* graph, path_handle_t ref_path,
                                                          const unordered_map<handle_t, vector<int64_t>>& handle_to_position,
                                                          int64_t max_deletion, int64_t context_steps,
                                                          const unordered_set<edge_t>& edge_blacklist) {
    vector<pair<int64_t, handle_t>> pos_handles;
    int64_t cur_pos = 0;
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            pos_handles.push_back(make_pair(cur_pos, handle));
            cur_pos += graph->get_length(handle);
        });

    vector<unordered_map<nid_t, int64_t>> id_to_min_pos_threads(get_thread_count());
    vector<unordered_map<nid_t, int64_t>> id_to_max_pos_threads(get_thread_count());
    
#pragma omp parallel for
    for (size_t i = 0; i < pos_handles.size(); ++i) {
        int64_t pos = pos_handles[i].first;
        handle_t handle = pos_handles[i].second;
        unordered_map<nid_t, int64_t>& id_to_min_pos = id_to_min_pos_threads[omp_get_thread_num()];
        unordered_map<nid_t, int64_t>& id_to_max_pos = id_to_max_pos_threads[omp_get_thread_num()];
        
        // scan a context of our current step, trying to avoid touching back
        // on the reference path
        // todo: can do better job of constraining with more step_on_path lookups
        unordered_set<nid_t> context;
        vector<handle_t> cur_handles = {handle};
        for (int64_t i = 0; i < context_steps && !cur_handles.empty(); ++i) {
            vector<handle_t> next_handles;
            for (auto& h : cur_handles) {
                nid_t cur_id = graph->get_id(h);
                if (!context.count(cur_id)) {
                    context.insert(cur_id);
                    graph->follow_edges(h, false, [&](handle_t n) {
                            if (!edge_blacklist.count(graph->edge_handle(h, n)) && !handle_to_position.count(n)) {
                                next_handles.push_back(n);
                            }
                        });
                    graph->follow_edges(h, true, [&](handle_t p) {
                            if (!edge_blacklist.count(graph->edge_handle(p, h)) && !handle_to_position.count(p)) {
                                next_handles.push_back(p);
                            }
                        });
                }
            }
            cur_handles = std::move(next_handles);
        }

        // assig everything in the context to the current position
        for (nid_t id : context) {
            auto it = id_to_min_pos.find(id);
            if (it == id_to_min_pos.end()) {
                id_to_min_pos[id] = pos;
                id_to_max_pos[id] = pos;
            } else {
                id_to_min_pos[id] = min(pos, it->second);
                id_to_max_pos[id] = max(pos, id_to_max_pos.at(id));
            }
        }                
    }

    for (size_t i = 1; i < id_to_max_pos_threads.size(); ++i) {
        for (const auto& id_max : id_to_max_pos_threads[i]) {
            id_to_max_pos_threads[0][id_max.first] = max(id_to_max_pos_threads[0][id_max.first], id_max.second);
        }
        id_to_max_pos_threads[i].clear();
        for (const auto& id_min : id_to_min_pos_threads[i]) {
            if (id_to_min_pos_threads[0].count(id_min.first)) {
                id_to_min_pos_threads[0][id_min.first] = min(id_to_min_pos_threads[0][id_min.first], id_min.second);
            } else {
                id_to_min_pos_threads[0][id_min.first] = id_min.second;
            }
        }
        id_to_min_pos_threads[i].clear();            
    }

    auto& id_to_min_pos = id_to_min_pos_threads[0];
    auto& id_to_max_pos = id_to_max_pos_threads[0];

    // scan every edge to find minimum distance according to positions found above
    multimap<int64_t, edge_t> length_to_edge;
    unordered_set<edge_t> edges_visited;
    vector<edge_t> neighbours;
    for (const auto& id_pos : id_to_min_pos) {
        handle_t handle = graph->get_handle(id_pos.first);
        graph->follow_edges(handle, false, [&] (handle_t next) {
                edge_t edge = graph->edge_handle(handle, next);
                if (id_to_min_pos.count(graph->get_id(next)) && !edges_visited.count(edge)) {
                    edges_visited.insert(edge);
                    neighbours.push_back(edge);
                }
            });
        graph->follow_edges(handle, true, [&] (handle_t next) {
                edge_t edge = graph->edge_handle(next, handle);
                if (id_to_min_pos.count(graph->get_id(next)) && !edges_visited.count(edge)) {
                    edges_visited.insert(edge);
                    neighbours.push_back(edge);
                }
            });
        for (const edge_t& edge : neighbours) {
            assert(graph->has_edge(edge));
            nid_t id1 = graph->get_id(edge.first);
            nid_t id2 = graph->get_id(edge.second);
            int64_t delta = max(abs(id_to_max_pos.at(id1) - id_to_min_pos.at(id2)), abs(id_to_max_pos.at(id2) - id_to_min_pos.at(id1)));
            
            if (delta > max_deletion) {
                length_to_edge.insert(make_pair(delta, edge));
            }
        }
        neighbours.clear();
    }

    return length_to_edge;
}

// walk from given edge back to reference path, returning the maximum and minimum distances found between
// the reference positions reached on its two sides (or -1s if a side reaches none)
static pair<int64_t, int64_t> get_deletion_span(const PathHandleGraph* graph, 
                                const unordered_map<handle_t, vector<int64_t>>& handle_to_position,
                                int64_t context_steps,
                                edge_t edge,
                                const unordered_set<edge_t>& edge_blacklist) {
    
    function<unordered_set<handle_t>(handle_t)> get_ref_context = [&] (handle_t handle) {
        unordered_set<handle_t> context;
        unordered_set<handle_t> ref_context;
        vector<handle_t> cur_handles = {handle};
        for (int64_t i = 0; i < context_steps && !cur_handles.empty(); ++i) {
            vector<handle_t> next_handles;
            for (auto& h : cur_handles) {
                if (!context.count(h)) {
                    context.insert(h);
                    if (i == 0) {
                        // keep search directional from origin
                        context.insert(graph->flip(h));
                    }
                    // stop search once we hit reference path
                    if (handle_to_position.count(h)) {
                        ref_context.insert(h);
                    } else {
                        graph->follow_edges(h, false, [&](handle_t n) {
                                if (!edge_blacklist.count(graph->edge_handle(h, n))) {
                                    next_handles.push_back(n);
                                }
                            });
                        if (i > 0) {
                            // keep search directional from origin
                            graph->follow_edges(h, true, [&](handle_t p) {
                                    if (!edge_blacklist.count(graph->edge_handle(p, h))) {
                                        next_handles.push_back(p);
                                    }
                                });
                        }
                    }
                }
            }
            cur_handles = std::move(next_handles);
        }
        return ref_context;
    };

    // search away from either end of the edge for the reference handles it reaches
    unordered_set<handle_t> left_context = get_ref_context(graph->flip(edge.first));
    unordered_set<handle_t> right_context = get_ref_context(edge.second);

    // a reference node that both searches reach is still connected to both sides without this edge, so
    // cutting the edge can't separate it from anything: leave it out. (otherwise every bubble along a long
    // deletion has an allele cut, though the other allele keeps the deletion connected)
    unordered_set<nid_t> left_ids, right_ids;
    for (handle_t handle : left_context) {
        left_ids.insert(graph->get_id(handle));
    }
    for (handle_t handle : right_context) {
        right_ids.insert(graph->get_id(handle));
    }
    vector<const vector<int64_t>*> left_positions, right_positions;
    bool revisited = false;
    for (handle_t handle : left_context) {
        if (!right_ids.count(graph->get_id(handle))) {
            left_positions.push_back(&handle_to_position.at(handle));
            revisited = revisited || left_positions.back()->size() > 1;
        }
    }
    for (handle_t handle : right_context) {
        if (!left_ids.count(graph->get_id(handle))) {
            right_positions.push_back(&handle_to_position.at(handle));
            revisited = revisited || right_positions.back()->size() > 1;
        }
    }
    if (left_positions.empty() || right_positions.empty()) {
        return make_pair(-1, -1);
    }

    // the minimum: the closest any visits on the two sides are
    vector<int64_t> left_sorted, right_sorted;
    for (const vector<int64_t>* positions : left_positions) {
        left_sorted.insert(left_sorted.end(), positions->begin(), positions->end());
    }
    for (const vector<int64_t>* positions : right_positions) {
        right_sorted.insert(right_sorted.end(), positions->begin(), positions->end());
    }
    sort(left_sorted.begin(), left_sorted.end());
    sort(right_sorted.begin(), right_sorted.end());
    int64_t min_delta = numeric_limits<int64_t>::max();
    for (size_t i = 0, j = 0; i < left_sorted.size() && j < right_sorted.size();) {
        min_delta = min(min_delta, abs(left_sorted[i] - right_sorted[j]));
        if (left_sorted[i] < right_sorted[j]) {
            ++i;
        } else {
            ++j;
        }
    }

    // compute the maximum deletion: the furthest apart any left and right reference handles are
    int64_t delta = 0;
    if (!revisited) {
        int64_t min_pos_left = numeric_limits<int64_t>::max(), max_pos_left = -1;
        int64_t min_pos_right = numeric_limits<int64_t>::max(), max_pos_right = -1;
        for (const vector<int64_t>* positions : left_positions) {
            min_pos_left = min(min_pos_left, positions->front());
            max_pos_left = max(max_pos_left, positions->front());
        }
        for (const vector<int64_t>* positions : right_positions) {
            min_pos_right = min(min_pos_right, positions->front());
            max_pos_right = max(max_pos_right, positions->front());
        }
        delta = max(abs(max_pos_left - min_pos_right), abs(max_pos_right - min_pos_left));
    } else {
        // a handle the reference visits more than once (it's collapsed) is taken at whichever visit is
        // closest to the other handle, since nothing says which copy a haplotype means
        for (const vector<int64_t>* left_pos : left_positions) {
            for (const vector<int64_t>* right_pos : right_positions) {
                int64_t closest = numeric_limits<int64_t>::max();
                for (int64_t l : *left_pos) {
                    for (int64_t r : *right_pos) {
                        closest = min(closest, abs(l - r));
                    }
                }
                delta = max(delta, closest);
            }
        }
    }
    
    return make_pair(delta, min_delta);
}

static int64_t edge_depth(PathHandleGraph* graph, edge_t edge, int64_t max_val) {
    // make absolutely sure it's canonical
    edge = graph->edge_handle(edge.first, edge.second);
    int64_t depth = 0;
    graph->for_each_step_on_handle(edge.first, [&](step_handle_t step) {
        path_handle_t path_handle = graph->get_path_handle_of_step(step);
        step_handle_t prev = graph->get_previous_step(step);
        if (prev != graph->path_front_end(path_handle) &&
            graph->edge_handle(graph->get_handle_of_step(prev), graph->get_handle_of_step(step)) == edge) {
            ++depth;
        } else {
            step_handle_t next = graph->get_next_step(step);
            if (next != graph->path_end(path_handle) &&
                graph->edge_handle(graph->get_handle_of_step(step), graph->get_handle_of_step(next)) == edge) {
                ++depth;
            }
        }
        return depth < max_val;
    });
    return depth;
}

void clip_deletion_edges(MutablePathMutableHandleGraph* graph, int64_t max_deletion,
                         int64_t context_steps, int64_t min_depth,
                         const vector<string>& ref_prefixes, int64_t min_fragment_len, bool verbose) {
    
    // load up the reference paths and their ids
    unordered_set<path_handle_t> ref_paths;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            for (const string& ref_prefix : ref_prefixes) {
                if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                    ref_paths.insert(path_handle);
                    break;
                }
            }
        });

    // all deletion edges for all paths
    unordered_set<edge_t> deletion_edges;
    // all reference edges + deletion edges
    unordered_set<edge_t> edge_blacklist;

    // the edges every reference path walks are never cut, and no search crosses them, so add them all
    // before searching from any one path (a path searched first would otherwise treat a later one as
    // off-reference sequence)
    unordered_set<edge_t> ref_edges;
    for (path_handle_t ref_path : ref_paths) {
        make_ref_index(graph, ref_path, ref_edges);
    }
    edge_blacklist = ref_edges;

    // order edges by node ids, so that ties don't depend on hashing, threads or graph implementation
    auto edge_key = [&](const edge_t& edge) {
        return make_tuple(graph->get_id(edge.first), graph->get_is_reverse(edge.first),
                          graph->get_id(edge.second), graph->get_is_reverse(edge.second));
    };
    
    for (path_handle_t ref_path : ref_paths) {

        // index the path to avoid path position interface dep
        unordered_map<handle_t, vector<int64_t>> handle_to_position = make_ref_index(graph, ref_path, edge_blacklist);

        // find set of deletion candidates sorted by length by walking out from the reference path
        if (verbose) {
            cerr << "[vg clip]: Searching for deletion candidates on " << graph->get_path_name(ref_path)
                 << " with " << context_steps << " context steps and " << get_thread_count() << " threads" << endl;
        }
        multimap<int64_t, edge_t> length_to_edge = find_deletion_candidate_edges(graph, ref_path, handle_to_position,
                                                                                 max_deletion, context_steps, edge_blacklist);
        // the candidates as (length, path depth, edge), leaving out reference edges, and edges with
        // depth >= min_depth if given
        vector<tuple<int64_t, int64_t, edge_t>> candidates;
        for (const auto& len_edge : length_to_edge) {
            if (!ref_edges.count(len_edge.second)) {
                int64_t depth = edge_depth(graph, len_edge.second, numeric_limits<int64_t>::max());
                if (min_depth <= 0 || depth < min_depth) {
                    candidates.emplace_back(len_edge.first, depth, len_edge.second);
                }
            }
        }
#ifdef debug
        cerr << "candidate edge count: " << length_to_edge.size() << " --> filtered count " << candidates.size() << endl;
#endif
        // longest first. the edges along one deletion tend to tie, and only one of them needs cutting, so
        // take the one the fewest paths cross: it's the deletion's own edge rather than an edge it shares
        // with haplotypes that stay local
        sort(candidates.begin(), candidates.end(), [&](const tuple<int64_t, int64_t, edge_t>& a,
                                                       const tuple<int64_t, int64_t, edge_t>& b) {
            if (get<0>(a) != get<0>(b)) {
                return get<0>(a) > get<0>(b);
            }
            if (get<1>(a) != get<1>(b)) {
                return get<1>(a) < get<1>(b);
            }
            return edge_key(get<2>(a)) < edge_key(get<2>(b));
        });

        if (verbose) {
            cerr << "[vg clip]: Found " << candidates.size() << " candidate deletion edges for " << graph->get_path_name(ref_path);
            if (!candidates.empty()) {
                cerr << " with sizes ranging from " << get<0>(candidates.back()) << " to " << get<0>(candidates.front());
            }
            cerr << endl;
        }

        // for every deletion candidate, walk *back* to the reference path and add it if it forms a deletion
        // we do this one-at-a-time to make sure that the edges are still deletions after the preceding edges were deleted.
        // and we do it twice: first only cutting edges that are nothing but a deletion (everything on one side
        // is far from everything on the other), then any edge that still forms one. so where a deletion runs
        // through an edge that local haplotypes also use, it's the deletion's own edge that goes
        for (bool clean_only : {true, false}) {
            int64_t candidate_i = 0;
            for (const auto& candidate : candidates) {
                const edge_t& edge = get<2>(candidate);
                assert(graph->has_edge(edge));
                if (!deletion_edges.count(edge)) {
                    pair<int64_t, int64_t> span = get_deletion_span(graph, handle_to_position, context_steps, edge, edge_blacklist);

                    if ((clean_only ? span.second : span.first) > max_deletion) {
                        deletion_edges.insert(edge);
                        edge_blacklist.insert(edge);
                        if (verbose) {
                            cerr << "[vg clip]: Found deletion edge for candidate " << candidate_i << " on " << graph->get_path_name(ref_path) << ": "
                                 << (graph->get_is_reverse(edge.first) ? "<" : ">") << graph->get_id(edge.first)
                                 << (graph->get_is_reverse(edge.second) ? "<" : ">") << graph->get_id(edge.second)
                                 << " with reference length " << span.first << endl;
                        }
                    }
                }
                ++candidate_i;
            }
        }
    }

    // just for logging
    unordered_map<string, size_t> clip_counts;

    if (verbose) {
        cerr << "[vg-clip]: Clipping " << deletion_edges.size() << " edges" << endl;
    }

    // delete the edges
    delete_nodes_and_chop_paths(graph, {}, deletion_edges, min_fragment_len, verbose ? &clip_counts : nullptr);

    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path " << kv.first << endl;
        }
        clip_counts.clear();
    }
}

void clip_stubs_generic(MutablePathMutableHandleGraph* graph,
                        function<void(function<void(handle_t, const Region*)>)> iterate_handles,
                        function<bool(handle_t)> handle_in_range,
                        const vector<string>& ref_prefixes,
                        int64_t min_fragment_len,
                        bool verbose) {
    
    unordered_set<nid_t> to_delete;

    // just for logging
    unordered_map<string, size_t> clip_counts;

    // frontier for recursing on stub neighbours
    unordered_map<handle_t, const Region*> stub_neighbours_1;

    // test if a node is "reference" using a name check
    function<bool(const string&)> check_prefixes = [&ref_prefixes] (const string& path_name) {
        for (const string& ref_prefix : ref_prefixes) {
            if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                return true;
            }
        }
        return false;
    };

    // test if a node is a stub.
    // we consider a node a stub if either (or both) sides have degree 0.
    function<bool(const handle_t& handle)> is_stub = [&to_delete, &graph] (const handle_t& handle) {
        size_t left_degree = 0;
        graph->follow_edges(handle, true, [&](handle_t left) {
            if (!to_delete.count(graph->get_id(left))) {
                ++left_degree;
                return false;
            }
            return true;
        });
        size_t right_degree = 1;
        if (left_degree > 0) {
            right_degree = 0;
            graph->follow_edges(handle, false, [&](handle_t right) {
                if (!to_delete.count(graph->get_id(right))) {
                    ++right_degree;
                    return false;
                }
                return true;
            });
        }
        return left_degree == 0 || right_degree == 0;
    };
    
    function<void(handle_t, const Region*)> visit_handle = [&](handle_t handle, const Region* region) {

        if (!to_delete.count(graph->get_id(handle)) && is_stub(handle)) {
            bool on_ref = false;
            graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle) {
                if (!ref_prefixes.empty() || region) {
                    // if we have a region, do exact comparison to it.
                    // otherwise, do a prefix check against ref_prefix
                    string path_name = graph->get_path_name(graph->get_path_handle_of_step(step_handle));
                    if ((region && region->seq == path_name) || (!region && check_prefixes(path_name))) {
                        on_ref = true;
                        return false;
                    }
                }
                return true;
            });
            if (!on_ref) {
                to_delete.insert(graph->get_id(handle));

                // remember the neighbours -- they can be new stubs!
                graph->follow_edges(handle, true, [&](handle_t prev) {
                    if (handle_in_range(prev) && !to_delete.count(graph->get_id(prev)) && graph->get_id(handle) != graph->get_id(prev)) {
                        stub_neighbours_1[prev] = region;
                    }
                });
                graph->follow_edges(handle, false, [&](handle_t next) {
                    if (handle_in_range(next) && !to_delete.count(graph->get_id(next)) && graph->get_id(handle) != graph->get_id(next)) {
                        stub_neighbours_1[next] = region;
                    }
                });
            }
        }
    };
       
    // first pass: find all the stubs in iterate_handles
    // and populate stub_neighbours_1
    iterate_handles(visit_handle);

    // keep doing the same thing on the neighbours until none left, using
    // handle_in_range to make sure we don't step out of bounds (in the case we're doing BED regions)           
    unordered_map<handle_t, const Region*> stub_neighbours_2;
    while (!stub_neighbours_1.empty()) {
        stub_neighbours_2.clear();
        swap(stub_neighbours_1, stub_neighbours_2);
        for (const auto&  neighbour_pair : stub_neighbours_2) {
            visit_handle(neighbour_pair.first, neighbour_pair.second);
        }
        stub_neighbours_2.clear();
    }
        
    if (verbose) {
        cerr << "[vg-clip]: Removing " << to_delete.size() << " nodes from (non-reference) stubs." << endl;
    }

    // cut out the nodes and chop up paths
    delete_nodes_and_chop_paths(graph, to_delete, {}, min_fragment_len, verbose ? &clip_counts : nullptr);
        
    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path " << kv.first << endl;
        }
        clip_counts.clear();
    }

}

void clip_stubs(MutablePathMutableHandleGraph* graph, const vector<string>& ref_prefixes, int64_t min_fragment_len, bool verbose) {
    
    function<void(function<void(handle_t, const Region*)>)> iterate_handles = [&] (function<void(handle_t, const Region*)> visit_handle) {
        graph->for_each_handle([&](handle_t handle) {
                visit_handle(handle, nullptr);
            });        
    };

    function<bool(handle_t)> handle_in_range = [](handle_t) {
        return true;
    };

    clip_stubs_generic(graph, iterate_handles, handle_in_range, ref_prefixes, min_fragment_len, verbose);
}

void clip_contained_stubs(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                          SnarlManager& snarl_manager, bool include_endpoints, int64_t min_fragment_len, bool verbose) {

    unordered_set<handle_t> all_handles;
    function<bool(handle_t)> handle_in_range = [&all_handles](handle_t handle) {
        return all_handles.count(handle);
    };
    
    function<void(function<void(handle_t, const Region*)>)> iterate_handles = [&] (function<void(handle_t, const Region*)> visit_handle) {
        
        visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                                                                                        int64_t start_pos, int64_t end_pos,
                                                                                        bool steps_reversed, const Region* containing_region) {
                                   
            pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
            for (id_t node_id : contents.first) {
                visit_handle(graph->get_handle(node_id), containing_region);
                all_handles.insert(graph->get_handle(node_id));
            }                                   
        });
    };

    // the edge depths are computed globally, without looking at regions.  as such, they need some notion of reference paths
    // so we shimmy a set in from the regions
    set<string> ref_path_set;
    for (const Region& region : regions) {
        ref_path_set.insert(region.seq);
    }
    vector<string> ref_paths_from_regions(ref_path_set.begin(), ref_path_set.end());

    clip_stubs_generic(graph, iterate_handles, handle_in_range, ref_paths_from_regions, min_fragment_len, verbose);
    
}

void stubbify_ref_paths(MutablePathMutableHandleGraph* graph, const vector<string>& ref_prefixes, int64_t min_fragment_len,
                        bool keep_connected, bool verbose) {
    // the selected reference paths, and the edges they walk. those are never cut, even where a path
    // comes back through one of its own ends, so a reference is never chopped
    vector<path_handle_t> ref_paths;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        for (const string& ref_prefix : ref_prefixes) {
            if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                if (!graph->is_empty(path_handle)) {
                    ref_paths.push_back(path_handle);
                }
                break;
            }
        }
    });
    unordered_set<edge_t> ref_edges;
    for (path_handle_t path_handle : ref_paths) {
        handle_t prev_handle;
        bool has_prev = false;
        graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            if (has_prev) {
                ref_edges.insert(graph->edge_handle(prev_handle, handle));
            }
            prev_handle = handle;
            has_prev = true;
        });
    }

    // every other edge on the outer side of a reference path end, once, as (end, next): end is the end
    // node's handle reading out of the path (the flipped first handle, or the last handle), and next is
    // the neighbour that follows it reading outward
    vector<pair<handle_t, handle_t>> outer_edges;
    vector<size_t> outer_edge_paths; // the index in ref_paths of the path each one was found on (just for logging)
    unordered_set<edge_t> outer_edge_set;
    unordered_set<nid_t> end_ids;
    for (size_t p = 0; p < ref_paths.size(); ++p) {
        handle_t first_out = graph->flip(graph->get_handle_of_step(graph->path_begin(ref_paths[p])));
        handle_t last_out = graph->get_handle_of_step(graph->path_back(ref_paths[p]));
        for (handle_t end : {first_out, last_out}) {
            end_ids.insert(graph->get_id(end));
            graph->follow_edges(end, false, [&](handle_t next) {
                edge_t edge = graph->edge_handle(end, next);
                if (!ref_edges.count(edge) && outer_edge_set.insert(edge).second) {
                    outer_edges.emplace_back(end, next);
                    outer_edge_paths.push_back(p);
                }
            });
        }
    }

    unordered_set<edge_t> edges_to_delete;
    // just for logging
    vector<bool> path_clipped(ref_paths.size(), false);
    size_t kept_piece_count = 0, kept_node_count = 0, kept_edge_count = 0, tipless_piece_count = 0;
    int64_t kept_bp = 0;

    if (!keep_connected) {
        edges_to_delete = outer_edge_set;
        for (size_t p : outer_edge_paths) {
            path_clipped[p] = true;
        }
    } else if (!outer_edges.empty()) {
        // an end can only be cut away from the reference by its own outer edges if those edges are all
        // that attach some piece of the graph. cutting the rest is enough to take the end off any cycle,
        // so find the pieces: flood from the reference paths without crossing outer edges
        auto is_outer = [&](const handle_t& from, const handle_t& to) {
            return (end_ids.count(graph->get_id(from)) || end_ids.count(graph->get_id(to))) &&
                outer_edge_set.count(graph->edge_handle(from, to));
        };
        // flood from the handles on the stack, pushing each node that visit() reports as new
        function<void(vector<handle_t>&, const function<bool(nid_t)>&)> flood = [&](vector<handle_t>& stack,
                                                                                   const function<bool(nid_t)>& visit) {
            while (!stack.empty()) {
                handle_t handle = stack.back();
                stack.pop_back();
                for (bool go_left : {false, true}) {
                    graph->follow_edges(handle, go_left, [&](handle_t next) {
                        if (!(go_left ? is_outer(next, handle) : is_outer(handle, next)) && visit(graph->get_id(next))) {
                            stack.push_back(graph->forward(next));
                        }
                    });
                }
            }
        };

        nid_t min_id = graph->min_node_id();
        vector<bool> anchored(graph->max_node_id() - min_id + 1, false);
        function<bool(nid_t)> visit_anchored = [&](nid_t node_id) {
            if (anchored[node_id - min_id]) {
                return false;
            }
            anchored[node_id - min_id] = true;
            return true;
        };
        vector<handle_t> stack;
        for (path_handle_t path_handle : ref_paths) {
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = graph->forward(graph->get_handle_of_step(step_handle));
                if (visit_anchored(graph->get_id(handle))) {
                    stack.push_back(handle);
                }
            });
        }
        flood(stack, visit_anchored);

        // the pieces that only outer edges attach, with the outer edges (as indexes) that attach them
        unordered_map<nid_t, size_t> piece_of;
        vector<vector<nid_t>> pieces;
        vector<vector<size_t>> piece_outer_edges;
        for (size_t i = 0; i < outer_edges.size(); ++i) {
            nid_t next_id = graph->get_id(outer_edges[i].second);
            if (anchored[next_id - min_id]) {
                // still attached to the reference without this edge, so it only closes a cycle
                edges_to_delete.insert(graph->edge_handle(outer_edges[i].first, outer_edges[i].second));
                path_clipped[outer_edge_paths[i]] = true;
                continue;
            }
            if (!piece_of.count(next_id)) {
                size_t piece = pieces.size();
                pieces.emplace_back();
                piece_outer_edges.emplace_back();
                function<bool(nid_t)> visit_piece = [&](nid_t node_id) {
                    if (piece_of.count(node_id)) {
                        return false;
                    }
                    piece_of[node_id] = piece;
                    pieces[piece].push_back(node_id);
                    return true;
                };
                visit_piece(next_id);
                stack.push_back(graph->get_handle(next_id));
                flood(stack, visit_piece);
            }
            piece_outer_edges[piece_of[next_id]].push_back(i);
        }

        // the path steps that cross an edge
        unordered_map<edge_t, int64_t> coverage_memo;
        auto coverage = [&](const edge_t& edge) {
            auto it = coverage_memo.find(edge);
            if (it == coverage_memo.end()) {
                it = coverage_memo.emplace(edge, edge_depth(graph, edge, numeric_limits<int64_t>::max())).first;
            }
            return it->second;
        };
        // the edges on one side of a node that aren't being cut
        auto live_edges = [&](const handle_t& handle, bool go_left) {
            vector<edge_t> edges;
            graph->follow_edges(handle, go_left, [&](handle_t other) {
                edge_t edge = go_left ? graph->edge_handle(other, handle) : graph->edge_handle(handle, other);
                if (!edges_to_delete.count(edge) && find(edges.begin(), edges.end(), edge) == edges.end()) {
                    edges.push_back(edge);
                }
            });
            return edges;
        };
        // order edges by node ids, so that ties don't depend on hashing or graph implementation
        auto edge_key = [&](const edge_t& edge) {
            return make_tuple(graph->get_id(edge.first), graph->get_is_reverse(edge.first),
                              graph->get_id(edge.second), graph->get_is_reverse(edge.second));
        };

        for (size_t piece = 0; piece < pieces.size(); ++piece) {
            // keep the piece on the end whose outer edges the paths cross most (the first such end on
            // ties). keeping its edges to two ends would close a cycle through both of them
            vector<pair<handle_t, int64_t>> end_coverage;
            for (size_t i : piece_outer_edges[piece]) {
                int64_t cov = coverage(graph->edge_handle(outer_edges[i].first, outer_edges[i].second));
                auto it = find_if(end_coverage.begin(), end_coverage.end(), [&](const pair<handle_t, int64_t>& end_cov) {
                    return end_cov.first == outer_edges[i].first;
                });
                if (it == end_coverage.end()) {
                    end_coverage.emplace_back(outer_edges[i].first, cov);
                } else {
                    it->second += cov;
                }
            }
            handle_t kept_end = max_element(end_coverage.begin(), end_coverage.end(),
                                            [](const pair<handle_t, int64_t>& a, const pair<handle_t, int64_t>& b) {
                                                return a.second < b.second;
                                            })->first;
            vector<size_t> kept;
            for (size_t i : piece_outer_edges[piece]) {
                if (outer_edges[i].first == kept_end) {
                    kept.push_back(i);
                } else {
                    edges_to_delete.insert(graph->edge_handle(outer_edges[i].first, outer_edges[i].second));
                    path_clipped[outer_edge_paths[i]] = true;
                }
            }

            bool has_tip = false;
            for (nid_t node_id : pieces[piece]) {
                handle_t handle = graph->get_handle(node_id);
                if (live_edges(handle, false).empty() || live_edges(handle, true).empty()) {
                    has_tip = true;
                    break;
                }
            }
            if (!has_tip) {
                // a piece with no tip of its own is all cycles, and hanging it off a reference end leaves the
                // snarl tree's root chain nowhere to end but a loop (vg haplotypes then rejects the graph).
                // one tip is enough, so cut every edge on one side of one node. to cut nothing off, take a
                // side that no edge of a spanning tree of the piece and its end uses (the tree prefers the
                // edges paths cross most), and of those the side whose edges paths cross least
                ++tipless_piece_count;
                unordered_map<nid_t, size_t> piece_index;
                for (size_t k = 0; k < pieces[piece].size(); ++k) {
                    piece_index[pieces[piece][k]] = k;
                }
                // vertex k is the piece's kth node, and the last one stands for the kept end, which is
                // all a live edge can leave the piece to
                size_t end_vertex = pieces[piece].size();
                vector<tuple<int64_t, edge_t, size_t, size_t>> piece_edges;
                unordered_set<edge_t> seen;
                for (nid_t node_id : pieces[piece]) {
                    handle_t handle = graph->get_handle(node_id);
                    for (bool go_left : {false, true}) {
                        graph->follow_edges(handle, go_left, [&](handle_t other) {
                            edge_t edge = go_left ? graph->edge_handle(other, handle) : graph->edge_handle(handle, other);
                            if (!edges_to_delete.count(edge) && seen.insert(edge).second) {
                                auto it = piece_index.find(graph->get_id(other));
                                piece_edges.emplace_back(coverage(edge), edge, piece_index.at(node_id),
                                                         it == piece_index.end() ? end_vertex : it->second);
                            }
                        });
                    }
                }
                sort(piece_edges.begin(), piece_edges.end(), [&](const tuple<int64_t, edge_t, size_t, size_t>& a,
                                                                 const tuple<int64_t, edge_t, size_t, size_t>& b) {
                    if (get<0>(a) != get<0>(b)) {
                        return get<0>(a) > get<0>(b);
                    }
                    return edge_key(get<1>(a)) < edge_key(get<1>(b));
                });
                structures::UnionFind components(end_vertex + 1, false);
                unordered_set<edge_t> tree;
                for (const auto& piece_edge : piece_edges) {
                    if (components.find_group(get<2>(piece_edge)) != components.find_group(get<3>(piece_edge))) {
                        components.union_groups(get<2>(piece_edge), get<3>(piece_edge));
                        tree.insert(get<1>(piece_edge));
                    }
                }
                // a leaf of the tree always has a side with no tree edge, so there is a side to cut
                vector<edge_t> best_side;
                int64_t best_cost = numeric_limits<int64_t>::max();
                nid_t best_id = 0;
                for (nid_t node_id : pieces[piece]) {
                    handle_t handle = graph->get_handle(node_id);
                    for (bool go_left : {true, false}) {
                        vector<edge_t> side = live_edges(handle, go_left);
                        int64_t cost = 0;
                        bool in_tree = false;
                        for (const edge_t& edge : side) {
                            cost += coverage(edge);
                            in_tree = in_tree || tree.count(edge);
                        }
                        if (!in_tree && (cost < best_cost || (cost == best_cost && node_id < best_id))) {
                            best_side = std::move(side);
                            best_cost = cost;
                            best_id = node_id;
                        }
                    }
                }
                assert(!best_side.empty());
                for (const edge_t& edge : best_side) {
                    edges_to_delete.insert(edge);
                }
                path_clipped[outer_edge_paths[kept.front()]] = true;
            }

            ++kept_piece_count;
            for (size_t i : kept) {
                if (!edges_to_delete.count(graph->edge_handle(outer_edges[i].first, outer_edges[i].second))) {
                    ++kept_edge_count;
                }
            }
            kept_node_count += pieces[piece].size();
            for (nid_t node_id : pieces[piece]) {
                kept_bp += graph->get_length(graph->get_handle(node_id));
            }
        }
    }

    // just for logging
    unordered_map<string, size_t> clip_counts;

    if (verbose) {
        if (keep_connected) {
            cerr << "[vg-clip]: Keeping " << kept_edge_count << " edges that attach " << kept_piece_count << " pieces ("
                 << kept_node_count << " nodes, " << kept_bp << " bp) to reference path ends, " << tipless_piece_count
                 << " of which had no tips and had one node side cut to make one" << endl;
        }
        cerr << "[vg-clip]: Clipping " << edges_to_delete.size() << " edges to stubbify "
             << count(path_clipped.begin(), path_clipped.end(), true) << " reference paths" << endl;
    }

    // delete the edges
    delete_nodes_and_chop_paths(graph, {}, edges_to_delete, min_fragment_len, verbose ? &clip_counts : nullptr);

    if (verbose) {
        for (const auto& kv : clip_counts) {
            cerr << "[vg-clip]: Ref path stubbification creating " << kv.second << " fragments from path " << kv.first << endl;
        }
        clip_counts.clear();
    }
}

}
