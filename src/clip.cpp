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

// the closest any visit in one sorted list of reference positions is to any in another
static int64_t closest_visits(const vector<int64_t>& a, const vector<int64_t>& b) {
    int64_t closest = numeric_limits<int64_t>::max();
    for (size_t i = 0, j = 0; i < a.size() && j < b.size();) {
        closest = min(closest, abs(a[i] - b[j]));
        if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return closest;
}

// the furthest apart any reference handle on the left and any on the right are (given as their lists of
// positions, which make_ref_index sorts), or -1 if a side has none: counting every handle (longest_shared),
// and counting only pairs where neither handle is shared (as flagged) with the other side (longest). a handle
// the reference visits more than once (it's collapsed) is taken at whichever visit is closest to the other
// handle, since nothing says which copy a haplotype means
static void longest_spans(const vector<pair<const vector<int64_t>*, bool>>& left_positions,
                          const vector<pair<const vector<int64_t>*, bool>>& right_positions,
                          int64_t& longest, int64_t& longest_shared) {
    longest = longest_shared = -1;
    if (left_positions.empty() || right_positions.empty()) {
        return;
    }
    bool revisited = false;
    for (const auto& positions : left_positions) {
        revisited = revisited || positions.first->size() > 1;
    }
    for (const auto& positions : right_positions) {
        revisited = revisited || positions.first->size() > 1;
    }
    if (!revisited) {
        // just the furthest apart ends
        for (bool count_shared : {true, false}) {
            int64_t min_pos_left = numeric_limits<int64_t>::max(), max_pos_left = -1;
            int64_t min_pos_right = numeric_limits<int64_t>::max(), max_pos_right = -1;
            for (const auto& positions : left_positions) {
                if (count_shared || !positions.second) {
                    min_pos_left = min(min_pos_left, positions.first->front());
                    max_pos_left = max(max_pos_left, positions.first->front());
                }
            }
            for (const auto& positions : right_positions) {
                if (count_shared || !positions.second) {
                    min_pos_right = min(min_pos_right, positions.first->front());
                    max_pos_right = max(max_pos_right, positions.first->front());
                }
            }
            if (max_pos_left >= 0 && max_pos_right >= 0) {
                (count_shared ? longest_shared : longest) = max(abs(max_pos_left - min_pos_right),
                                                                 abs(max_pos_right - min_pos_left));
            }
        }
        return;
    }
    for (const auto& left_pos : left_positions) {
        for (const auto& right_pos : right_positions) {
            int64_t delta = closest_visits(*left_pos.first, *right_pos.first);
            longest_shared = max(longest_shared, delta);
            if (!left_pos.second && !right_pos.second) {
                longest = max(longest, delta);
            }
        }
    }
}

// how far apart an edge's two sides are on the reference (all -1 if a side reaches none of it)
struct DeletionSpan {
    // leaving out reference nodes that both sides reach: the closest and furthest apart the two sides are
    int64_t shortest;
    int64_t longest;
    // and counting them
    int64_t longest_shared;
};

// walk from given edge back to reference path, measuring how far apart the reference positions reached on
// its two sides are. the searches never cross the edge itself, but can come back into its nodes the other
// way around, as where a haplotype turns around in a hairpin there
static DeletionSpan get_deletion_span(const PathHandleGraph* graph,
                                      const unordered_map<handle_t, vector<int64_t>>& handle_to_position,
                                      int64_t context_steps,
                                      edge_t edge,
                                      const unordered_set<edge_t>& edge_blacklist) {

    edge = graph->edge_handle(edge.first, edge.second);
    function<unordered_set<handle_t>(handle_t)> get_ref_context = [&] (handle_t handle) {
        unordered_set<handle_t> context;
        unordered_set<handle_t> ref_context;
        vector<handle_t> cur_handles = {handle};
        for (int64_t i = 0; i < context_steps && !cur_handles.empty(); ++i) {
            vector<handle_t> next_handles;
            for (auto& h : cur_handles) {
                if (!context.count(h)) {
                    context.insert(h);
                    // stop search once we hit reference path
                    if (handle_to_position.count(h)) {
                        ref_context.insert(h);
                    } else {
                        graph->follow_edges(h, false, [&](handle_t n) {
                                edge_t next_edge = graph->edge_handle(h, n);
                                if (next_edge != edge && !edge_blacklist.count(next_edge)) {
                                    next_handles.push_back(n);
                                }
                            });
                        if (i > 0) {
                            // keep search directional from origin
                            graph->follow_edges(h, true, [&](handle_t p) {
                                    edge_t prev_edge = graph->edge_handle(p, h);
                                    if (prev_edge != edge && !edge_blacklist.count(prev_edge)) {
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
    // cutting the edge can't separate it from anything. (counting it, every bubble along a long deletion
    // looks like a deletion edge, though the other allele keeps the deletion connected)
    unordered_set<nid_t> left_ids, right_ids;
    for (handle_t handle : left_context) {
        left_ids.insert(graph->get_id(handle));
    }
    for (handle_t handle : right_context) {
        right_ids.insert(graph->get_id(handle));
    }
    // each handle's positions, and whether the other side reaches its node too
    vector<pair<const vector<int64_t>*, bool>> left_positions, right_positions;
    vector<int64_t> left_sorted, right_sorted;
    for (handle_t handle : left_context) {
        left_positions.emplace_back(&handle_to_position.at(handle), right_ids.count(graph->get_id(handle)));
        if (!left_positions.back().second) {
            left_sorted.insert(left_sorted.end(), left_positions.back().first->begin(), left_positions.back().first->end());
        }
    }
    for (handle_t handle : right_context) {
        right_positions.emplace_back(&handle_to_position.at(handle), left_ids.count(graph->get_id(handle)));
        if (!right_positions.back().second) {
            right_sorted.insert(right_sorted.end(), right_positions.back().first->begin(), right_positions.back().first->end());
        }
    }

    DeletionSpan span;
    longest_spans(left_positions, right_positions, span.longest, span.longest_shared);
    span.shortest = -1;
    if (span.longest >= 0) {
        sort(left_sorted.begin(), left_sorted.end());
        sort(right_sorted.begin(), right_sorted.end());
        span.shortest = closest_visits(left_sorted, right_sorted);
    }
    return span;
}

// order edges by node ids, so that ties don't depend on hashing, threads or graph implementation
static tuple<nid_t, bool, nid_t, bool> edge_key(const HandleGraph* graph, const edge_t& edge) {
    return make_tuple(graph->get_id(edge.first), graph->get_is_reverse(edge.first),
                      graph->get_id(edge.second), graph->get_is_reverse(edge.second));
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
    
    // load up the reference paths and their ids, in name order so the results don't depend on the graph implementation
    vector<path_handle_t> ref_paths;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            for (const string& ref_prefix : ref_prefixes) {
                if (path_name.compare(0, ref_prefix.length(), ref_prefix) == 0) {
                    ref_paths.push_back(path_handle);
                    break;
                }
            }
        });
    sort(ref_paths.begin(), ref_paths.end(), [&](const path_handle_t& a, const path_handle_t& b) {
        return graph->get_path_name(a) < graph->get_path_name(b);
    });

    // all deletion edges for all paths
    unordered_set<edge_t> deletion_edges;

    // the edges any reference path walks are never cut. (while searching from one path, the others count
    // as off-reference sequence, so a jump through one of them can still be found. but it can't be removed:
    // where two selected paths share a component, an edge beside such a jump can be cut to no effect)
    unordered_set<edge_t> ref_edges;
    for (path_handle_t ref_path : ref_paths) {
        handle_t prev_handle;
        bool has_prev = false;
        graph->for_each_step_in_path(ref_path, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            if (has_prev) {
                ref_edges.insert(graph->edge_handle(prev_handle, handle));
            }
            prev_handle = handle;
            has_prev = true;
        });
    }
    
    for (path_handle_t ref_path : ref_paths) {

        // this path's reference edges + deletion edges, which no search crosses
        unordered_set<edge_t> edge_blacklist = deletion_edges;

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
        // the edges the fewest paths cross first: of the edges along one deletion, only one needs cutting,
        // and the one fewest haplotypes use is the deletion's own edge rather than one it shares with
        // haplotypes that stay local. then longest first
        sort(candidates.begin(), candidates.end(), [&](const tuple<int64_t, int64_t, edge_t>& a,
                                                       const tuple<int64_t, int64_t, edge_t>& b) {
            if (get<1>(a) != get<1>(b)) {
                return get<1>(a) < get<1>(b);
            }
            if (get<0>(a) != get<0>(b)) {
                return get<0>(a) > get<0>(b);
            }
            return edge_key(graph, get<2>(a)) < edge_key(graph, get<2>(b));
        });

        if (verbose) {
            cerr << "[vg clip]: Found " << candidates.size() << " candidate deletion edges for " << graph->get_path_name(ref_path);
            if (!candidates.empty()) {
                auto length_range = minmax_element(candidates.begin(), candidates.end(),
                                                   [](const tuple<int64_t, int64_t, edge_t>& a, const tuple<int64_t, int64_t, edge_t>& b) {
                                                       return get<0>(a) < get<0>(b);
                                                   });
                cerr << " with sizes ranging from " << get<0>(*length_range.first) << " to " << get<0>(*length_range.second);
            }
            cerr << endl;
        }

        auto edge_string = [&](const edge_t& edge) {
            return string(graph->get_is_reverse(edge.first) ? "<" : ">") + to_string(graph->get_id(edge.first)) +
                (graph->get_is_reverse(edge.second) ? "<" : ">") + to_string(graph->get_id(edge.second));
        };

        // each candidate's span, kept until the cuts change (with the number of changes when it was measured)
        size_t cut_changes = 0;
        vector<pair<size_t, DeletionSpan>> spans(candidates.size(), make_pair(numeric_limits<size_t>::max(), DeletionSpan()));
        auto measure = [&](size_t i) -> const DeletionSpan& {
            if (spans[i].first != cut_changes) {
                spans[i] = make_pair(cut_changes, get_deletion_span(graph, handle_to_position, context_steps,
                                                                    get<2>(candidates[i]), edge_blacklist));
            }
            return spans[i].second;
        };

        // for every deletion candidate, walk *back* to the reference path and add it if it forms a deletion
        // we do this one-at-a-time to make sure that the edges are still deletions after the preceding edges were deleted.
        // a sweep cuts the candidates whose span (as measured by the given function) is longer than max_deletion
        vector<size_t> cuts;
        auto sweep = [&](const function<int64_t(const DeletionSpan&)>& span_length) {
            for (size_t i = 0; i < candidates.size(); ++i) {
                const edge_t& edge = get<2>(candidates[i]);
                assert(graph->has_edge(edge));
                if (!deletion_edges.count(edge)) {
                    int64_t length = span_length(measure(i));
                    if (length > max_deletion) {
                        deletion_edges.insert(edge);
                        edge_blacklist.insert(edge);
                        ++cut_changes;
                        cuts.push_back(i);
                        if (verbose) {
                            cerr << "[vg clip]: Found deletion edge for candidate " << i << " on " << graph->get_path_name(ref_path) << ": "
                                 << edge_string(edge) << " with reference length " << length << endl;
                        }
                    }
                }
            }
        };
        // first, only edges that are nothing but a deletion: everything on one side is far from everything on
        // the other. so where a deletion runs through an edge that local haplotypes also use, it's the
        // deletion's own edge that goes
        sweep([](const DeletionSpan& span) { return span.shortest; });
        // then any edge that still forms one, first not counting the reference nodes both sides reach, then
        // counting them (as when a haplotype loops back to the same place at two far loci)
        sweep([](const DeletionSpan& span) { return span.longest; });
        sweep([](const DeletionSpan& span) { return span.longest_shared; });

        // cuts made later can leave earlier ones unnecessary: put back each one, latest first, that no longer
        // forms a deletion with the others in place
        vector<size_t> swept_cuts = std::move(cuts);
        vector<size_t> put_back;
        for (auto it = swept_cuts.rbegin(); it != swept_cuts.rend(); ++it) {
            const edge_t& edge = get<2>(candidates[*it]);
            edge_blacklist.erase(edge);
            if (get_deletion_span(graph, handle_to_position, context_steps, edge, edge_blacklist).longest_shared <= max_deletion) {
                deletion_edges.erase(edge);
                ++cut_changes;
                put_back.push_back(*it);
                if (verbose) {
                    cerr << "[vg clip]: Putting back deletion edge on " << graph->get_path_name(ref_path) << ": "
                         << edge_string(edge) << ", which later cuts made unnecessary" << endl;
                }
            } else {
                edge_blacklist.insert(edge);
            }
        }
        // and since putting edges back can join a deletion up again, one last sweep
        cuts.clear();
        sweep([](const DeletionSpan& span) { return span.longest_shared; });

        if (!cuts.empty()) {
            // the search overestimates where it can go, so with edges put back it can find jumps that no
            // haplotype makes, and cut more than was put back. the cuts from before putting any back left no
            // deletion (cutting more only shortens spans), so go back to those if they break fewer haplotypes
            auto cost = [&](const vector<size_t>& cut_set) {
                int64_t depth = 0;
                for (size_t i : cut_set) {
                    depth += get<1>(candidates[i]);
                }
                return make_pair(depth, cut_set.size());
            };
            unordered_set<size_t> put_back_set(put_back.begin(), put_back.end());
            vector<size_t> final_cuts = cuts;
            for (size_t i : swept_cuts) {
                if (!put_back_set.count(i)) {
                    final_cuts.push_back(i);
                }
            }
            if (cost(swept_cuts) < cost(final_cuts)) {
                for (size_t i : cuts) {
                    deletion_edges.erase(get<2>(candidates[i]));
                    edge_blacklist.erase(get<2>(candidates[i]));
                }
                for (size_t i : put_back) {
                    deletion_edges.insert(get<2>(candidates[i]));
                    edge_blacklist.insert(get<2>(candidates[i]));
                }
                if (verbose) {
                    cerr << "[vg clip]: Undoing " << put_back.size() << " put back and " << cuts.size()
                         << " later deletion edges on " << graph->get_path_name(ref_path)
                         << ", since the cuts from before putting any back break fewer paths" << endl;
                }
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
    // in name order, so that ties between paths don't depend on the graph implementation
    sort(ref_paths.begin(), ref_paths.end(), [&](const path_handle_t& a, const path_handle_t& b) {
        return graph->get_path_name(a) < graph->get_path_name(b);
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
    unordered_map<edge_t, size_t> outer_edge_index;
    vector<vector<size_t>> outer_edge_paths; // the indexes in ref_paths of the paths it's at the end of (just for logging)
    unordered_set<nid_t> end_ids;
    for (size_t p = 0; p < ref_paths.size(); ++p) {
        handle_t first_out = graph->flip(graph->get_handle_of_step(graph->path_begin(ref_paths[p])));
        handle_t last_out = graph->get_handle_of_step(graph->path_back(ref_paths[p]));
        for (handle_t end : {first_out, last_out}) {
            end_ids.insert(graph->get_id(end));
            graph->follow_edges(end, false, [&](handle_t next) {
                edge_t edge = graph->edge_handle(end, next);
                if (!ref_edges.count(edge)) {
                    auto it = outer_edge_index.find(edge);
                    if (it == outer_edge_index.end()) {
                        outer_edge_index[edge] = outer_edges.size();
                        outer_edges.emplace_back(end, next);
                        outer_edge_paths.push_back({p});
                    } else if (outer_edge_paths[it->second].back() != p) {
                        outer_edge_paths[it->second].push_back(p);
                    }
                }
            });
        }
    }

    unordered_set<edge_t> edges_to_delete;
    // just for logging
    vector<bool> path_clipped(ref_paths.size(), false);
    auto clip_outer_edge = [&](size_t i) {
        edges_to_delete.insert(graph->edge_handle(outer_edges[i].first, outer_edges[i].second));
        for (size_t p : outer_edge_paths[i]) {
            path_clipped[p] = true;
        }
    };
    size_t kept_piece_count = 0, kept_node_count = 0, kept_edge_count = 0, looped_end_count = 0, kept_ends_count = 0;
    int64_t kept_bp = 0;

    if (!keep_connected) {
        for (size_t i = 0; i < outer_edges.size(); ++i) {
            clip_outer_edge(i);
        }
    } else if (!outer_edges.empty()) {
        // an end can only be cut away from the reference by its own outer edges if those edges are all
        // that attach some piece of the graph. cutting the rest is enough to take the end off any cycle,
        // so find the pieces: flood from the reference paths without crossing outer edges
        auto is_outer = [&](const handle_t& from, const handle_t& to) {
            return (end_ids.count(graph->get_id(from)) || end_ids.count(graph->get_id(to))) &&
                outer_edge_index.count(graph->edge_handle(from, to));
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
                clip_outer_edge(i);
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
        // the ends pieces are kept on, in the order first used, with the pieces and their kept outer edges
        vector<handle_t> kept_ends;
        vector<vector<size_t>> end_pieces;
        vector<vector<size_t>> end_kept;
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
            auto end_it = find(kept_ends.begin(), kept_ends.end(), kept_end);
            if (end_it == kept_ends.end()) {
                end_it = kept_ends.insert(kept_ends.end(), kept_end);
                end_pieces.emplace_back();
                end_kept.emplace_back();
            }
            size_t end_number = end_it - kept_ends.begin();
            end_pieces[end_number].push_back(piece);
            for (size_t i : piece_outer_edges[piece]) {
                if (outer_edges[i].first == kept_end) {
                    end_kept[end_number].push_back(i);
                } else {
                    clip_outer_edge(i);
                }
            }
        }
        kept_ends_count = kept_ends.size();

        for (size_t end_number = 0; end_number < kept_ends.size(); ++end_number) {
            // vg haplotypes needs the snarl tree's root chain to run from the reference end out to a tip. it can't
            // if what hangs off the end, or a dead-end branch of it, is all cycles: the chain can end in one of
            // those instead. so every dead end gets a tip. look at the node sides of all the pieces on the end:
            // each node joins its two sides, and each edge joins two sides. split that at its bridges, and the
            // parts at the leaves are the dead ends (all of it if there are no bridges). a dead end that is just
            // one side with no edges is a tip; any other has cycles, and gets a tip by cutting every edge on one
            // of its sides. to cut nothing off, that's a side whose only edge in a spanning tree is its node (the
            // tree prefers the edges paths cross most; a leaf of it in the dead end is such a side), and of those
            // the side whose edges paths cross least, then the one furthest into the dead end (cutting nearer the
            // way in can leave the rest as another dead end with cycles). a cut can still leave one, so repeat
            // until there are none. (the pieces on an end go together so that one with a tip counts for all of
            // them, as it would once a cut splits a piece into two on the same end)
            const vector<size_t>& kept = end_kept[end_number];
            vector<nid_t> end_nodes;
            for (size_t piece : end_pieces[end_number]) {
                end_nodes.insert(end_nodes.end(), pieces[piece].begin(), pieces[piece].end());
            }
            unordered_map<nid_t, size_t> piece_index;
            for (size_t k = 0; k < end_nodes.size(); ++k) {
                piece_index[end_nodes[k]] = k;
            }
            // vertex 2k is the left side of the kth node and 2k + 1 its right side, and the last vertex stands for
            // the kept end's outer side, which is all a live edge can leave the pieces to
            size_t end_vertex = 2 * end_nodes.size();
            function<size_t(const handle_t&, bool)> side_of = [&](const handle_t& handle, bool right) {
                auto it = piece_index.find(graph->get_id(handle));
                return it == piece_index.end() ? end_vertex : 2 * it->second + (right != graph->get_is_reverse(handle) ? 1 : 0);
            };
            bool made_tip = false;
            while (true) {
                // the live edges, as (coverage, edge, side, side), well covered first
                vector<tuple<int64_t, edge_t, size_t, size_t>> piece_edges;
                unordered_set<edge_t> seen;
                for (nid_t node_id : end_nodes) {
                    handle_t handle = graph->get_handle(node_id);
                    for (bool go_left : {false, true}) {
                        graph->follow_edges(handle, go_left, [&](handle_t other) {
                            edge_t edge = go_left ? graph->edge_handle(other, handle) : graph->edge_handle(handle, other);
                            if (!edges_to_delete.count(edge) && seen.insert(edge).second) {
                                if (go_left) {
                                    piece_edges.emplace_back(coverage(edge), edge, side_of(other, true), side_of(handle, false));
                                } else {
                                    piece_edges.emplace_back(coverage(edge), edge, side_of(handle, true), side_of(other, false));
                                }
                            }
                        });
                    }
                }
                sort(piece_edges.begin(), piece_edges.end(), [&](const tuple<int64_t, edge_t, size_t, size_t>& a,
                                                                 const tuple<int64_t, edge_t, size_t, size_t>& b) {
                    if (get<0>(a) != get<0>(b)) {
                        return get<0>(a) > get<0>(b);
                    }
                    return edge_key(graph, get<1>(a)) < edge_key(graph, get<1>(b));
                });
                // the joins: the edges, then the nodes, numbered after them
                size_t join_count = piece_edges.size() + end_nodes.size();
                auto join_sides = [&](size_t j) {
                    return j < piece_edges.size() ? make_pair(get<2>(piece_edges[j]), get<3>(piece_edges[j])) :
                        make_pair(2 * (j - piece_edges.size()), 2 * (j - piece_edges.size()) + 1);
                };

                // the spanning tree, which takes every node first
                structures::UnionFind tree_components(end_vertex + 1, false);
                for (size_t k = 0; k < end_nodes.size(); ++k) {
                    tree_components.union_groups(2 * k, 2 * k + 1);
                }
                vector<bool> in_tree(piece_edges.size(), false);
                for (size_t e = 0; e < piece_edges.size(); ++e) {
                    if (tree_components.find_group(get<2>(piece_edges[e])) != tree_components.find_group(get<3>(piece_edges[e]))) {
                        tree_components.union_groups(get<2>(piece_edges[e]), get<3>(piece_edges[e]));
                        in_tree[e] = true;
                    }
                }

                // find the bridges, with a depth first search from the end (self loops can't be bridges)
                vector<vector<pair<size_t, size_t>>> adjacent(end_vertex + 1);
                for (size_t j = 0; j < join_count; ++j) {
                    pair<size_t, size_t> sides = join_sides(j);
                    if (sides.first != sides.second) {
                        adjacent[sides.first].emplace_back(sides.second, j);
                        adjacent[sides.second].emplace_back(sides.first, j);
                    }
                }
                vector<bool> is_bridge(join_count, false);
                size_t bridge_count = 0;
                vector<size_t> order(end_vertex + 1, 0), low(end_vertex + 1, 0);
                size_t visited = 0;
                // (vertex, the join it was reached by, the next adjacency to look at)
                vector<tuple<size_t, size_t, size_t>> search = {make_tuple(end_vertex, numeric_limits<size_t>::max(), 0)};
                order[end_vertex] = low[end_vertex] = ++visited;
                while (!search.empty()) {
                    size_t v = get<0>(search.back());
                    if (get<2>(search.back()) < adjacent[v].size()) {
                        pair<size_t, size_t> next = adjacent[v][get<2>(search.back())++];
                        if (next.second == get<1>(search.back())) {
                            continue;
                        }
                        if (order[next.first]) {
                            low[v] = min(low[v], order[next.first]);
                        } else {
                            order[next.first] = low[next.first] = ++visited;
                            search.emplace_back(next.first, next.second, 0);
                        }
                    } else {
                        size_t j = get<1>(search.back());
                        search.pop_back();
                        if (!search.empty()) {
                            size_t u = get<0>(search.back());
                            low[u] = min(low[u], low[v]);
                            if (low[v] > order[u]) {
                                is_bridge[j] = true;
                                ++bridge_count;
                            }
                        }
                    }
                }

                // the parts the bridges join, how many bridges each has, and whether each has an edge (if
                // not, it's a lone side with no edges: a tip)
                structures::UnionFind parts(end_vertex + 1, false);
                for (size_t j = 0; j < join_count; ++j) {
                    if (!is_bridge[j]) {
                        pair<size_t, size_t> sides = join_sides(j);
                        parts.union_groups(sides.first, sides.second);
                    }
                }
                vector<size_t> part_bridges(end_vertex + 1, 0);
                vector<bool> part_has_edge(end_vertex + 1, false);
                for (size_t j = 0; j < join_count; ++j) {
                    pair<size_t, size_t> sides = join_sides(j);
                    if (is_bridge[j]) {
                        ++part_bridges[parts.find_group(sides.first)];
                        ++part_bridges[parts.find_group(sides.second)];
                    } else if (j < piece_edges.size()) {
                        part_has_edge[parts.find_group(sides.first)] = true;
                    }
                }
                size_t end_part = parts.find_group(end_vertex);
                auto is_looped = [&](size_t part) {
                    return part_has_edge[part] && (part == end_part ? bridge_count == 0 : part_bridges[part] == 1);
                };

                // how far into its dead end each side is, in joins from the side the way in comes to
                vector<size_t> depth(end_vertex + 1, numeric_limits<size_t>::max());
                vector<size_t> queue;
                if (is_looped(end_part)) {
                    depth[end_vertex] = 0;
                    queue.push_back(end_vertex);
                }
                for (size_t j = 0; j < join_count; ++j) {
                    if (is_bridge[j]) {
                        pair<size_t, size_t> sides = join_sides(j);
                        for (size_t side : {sides.first, sides.second}) {
                            if (is_looped(parts.find_group(side)) && depth[side] != 0) {
                                depth[side] = 0;
                                queue.push_back(side);
                            }
                        }
                    }
                }
                for (size_t q = 0; q < queue.size(); ++q) {
                    for (const pair<size_t, size_t>& next : adjacent[queue[q]]) {
                        if (!is_bridge[next.second] && depth[next.first] == numeric_limits<size_t>::max()) {
                            depth[next.first] = depth[queue[q]] + 1;
                            queue.push_back(next.first);
                        }
                    }
                }

                // the side to cut in each dead end with cycles
                map<size_t, tuple<int64_t, size_t, vector<edge_t>>> cuts;
                unordered_set<size_t> looped_parts;
                vector<int64_t> side_cost(end_vertex, 0);
                vector<bool> side_in_tree(end_vertex, false);
                vector<vector<edge_t>> side_edges(end_vertex);
                for (size_t e = 0; e < piece_edges.size(); ++e) {
                    for (size_t side : {get<2>(piece_edges[e]), get<3>(piece_edges[e])}) {
                        if (side != end_vertex && (side_edges[side].empty() || side_edges[side].back() != get<1>(piece_edges[e]))) {
                            side_edges[side].push_back(get<1>(piece_edges[e]));
                            side_cost[side] += get<0>(piece_edges[e]);
                            side_in_tree[side] = side_in_tree[side] || in_tree[e];
                        }
                    }
                }
                // cheapest first, then furthest in, then lowest node id
                auto side_order = [&](size_t side) {
                    return make_tuple(side_cost[side], -(int64_t)depth[side], end_nodes[side / 2]);
                };
                for (size_t side = 0; side < end_vertex; ++side) {
                    size_t part = parts.find_group(side);
                    if (!is_looped(part)) {
                        continue;
                    }
                    looped_parts.insert(part);
                    auto it = cuts.find(part);
                    if (!side_in_tree[side] && !side_edges[side].empty() &&
                        (it == cuts.end() || side_order(side) < side_order(get<1>(it->second)))) {
                        cuts[part] = make_tuple(side_cost[side], side, side_edges[side]);
                    }
                }
                assert(cuts.size() == looped_parts.size());
                if (cuts.empty()) {
                    break;
                }
                for (const auto& part_cut : cuts) {
                    for (const edge_t& edge : get<2>(part_cut.second)) {
                        edges_to_delete.insert(edge);
                    }
                }
                made_tip = true;
            }
            if (made_tip) {
                ++looped_end_count;
                for (size_t p : outer_edge_paths[kept.front()]) {
                    path_clipped[p] = true;
                }
            }

            kept_piece_count += end_pieces[end_number].size();
            for (size_t i : kept) {
                if (!edges_to_delete.count(graph->edge_handle(outer_edges[i].first, outer_edges[i].second))) {
                    ++kept_edge_count;
                }
            }
            kept_node_count += end_nodes.size();
            for (nid_t node_id : end_nodes) {
                kept_bp += graph->get_length(graph->get_handle(node_id));
            }
        }
    }

    // just for logging
    unordered_map<string, size_t> clip_counts;

    if (verbose) {
        if (keep_connected) {
            cerr << "[vg-clip]: Keeping " << kept_edge_count << " edges that attach " << kept_piece_count << " pieces ("
                 << kept_node_count << " nodes, " << kept_bp << " bp) to " << kept_ends_count
                 << " reference path ends, on " << looped_end_count
                 << " of which dead ends with no tips had a node side cut in each to make one" << endl;
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
