#include "clip.hpp"
#include "traversal_finder.hpp"
#include <unordered_map>
#include <IntervalTree.h>
#include <structures/rank_pairing_heap.hpp>
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
    PathTraversalFinder trav_finder(*graph, snarl_manager, path_names);
    
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
                                 unordered_map<string, size_t>* fragments_per_path) {
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
                    if (step_deleted || cur_step_offset - start_step_offset >= min_fragment_len) {
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
            graph->destroy_path(path_handle);
        }
    }

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
                             size_t max_nodes, size_t max_edges,
                             size_t max_nodes_shallow, size_t max_edges_shallow,
                             double max_avg_degree, double max_reflen_prop, size_t max_reflen,
                             double& out_avg_degree) {

    out_avg_degree = -1.;
    
    // if our snarl is to big vs the reference path, we do not process it
    double ref_prop = (double)ref_interval_length / (double)graph->get_path_length(path_handle);
    if (ref_prop > max_reflen_prop || ref_interval_length > max_reflen) {
#ifdef debug
        cerr << "skipping snarl " << pb2json(*snarl) << " with interval length " << ref_interval_length
             << " because its ref_prop of " << region.seq << " is " << ref_prop << " which is greater than " << max_reflen_prop
             << " or its ref length " << ref_interval_length << " is greater than " << max_reflen << endl;
#endif
        return false;
    }

    // check the stats
    bool filter_on = max_nodes > 0 || max_edges > 0 || max_nodes_shallow > 0 || max_edges_shallow > 0 || max_avg_degree > 0.;
    bool complex_nodes = contents.first.size() > max_nodes;
    bool complex_edges = contents.second.size() > max_edges;
    bool complex_nodes_shallow = contents_shallow.first.size() > max_nodes_shallow;
    bool complex_edges_shallow = contents_shallow.second.size() > max_edges_shallow;
    size_t total_degree = 0;
    for (id_t node_id : contents.first) {
        handle_t handle = graph->get_handle(node_id);
        total_degree += graph->get_degree(handle, true) + graph->get_degree(handle, false);
    }
    // degree averaged over node sides to be a bit more intuitive, hence 2X in denominator:
    out_avg_degree = (double)total_degree / (2. *(double)contents.first.size());
    bool complex_degree = out_avg_degree > max_avg_degree;        
    
    return !filter_on || (complex_nodes && complex_edges && complex_nodes_shallow && complex_edges_shallow && complex_degree);
}

void clip_contained_snarls(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions, 
                           SnarlManager& snarl_manager, bool include_endpoints, int64_t min_fragment_len,
                           size_t max_nodes, size_t max_edges, size_t max_nodes_shallow, size_t max_edges_shallow,
                           double max_avg_degree, double max_reflen_prop, size_t max_reflen,
                           bool out_bed, bool verbose) {

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
            bool deletion_on_whitellist = whitelist.size() <= 2;
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
            if (snarl_is_complex(pp_graph, snarl, contents, contents_shallow, ref_interval_length, *containing_region, path_handle, max_nodes, max_edges,
                                 max_nodes_shallow, max_edges_shallow, max_avg_degree, max_reflen_prop, max_reflen, avg_degree)) {
                if (out_bed) {
                    string snarl_name = (snarl->start().backward() ? "<" : ">") + std::to_string(snarl->start().node_id()) +
                        (snarl->end().backward() ? "<" : ">") + std::to_string(snarl->end().node_id());
                    cout << containing_region->seq << "\t" << start_pos << "\t" << (end_pos + 1) << "\t" << snarl_name << "\t"
                         << contents.first.size() << "\t" << contents.second.size() << "\t"
                         << contents_shallow.first.size() << "\t" << contents_shallow.second.size() << "\t"
                         << avg_degree << "\n";
                } else {
                    for (id_t node_id : contents.first) {
                        if (!whitelist.count(node_id)) {
                            nodes_to_delete.insert(node_id);
                            ++clip_counts[containing_region->seq];
                        }
                    }
                    // since we're deleting all alt alleles, the only edge that could be left is a snarl-spanning deletion
                    if (!deletion_on_whitellist) {
                        edge_t deletion_edge = graph->edge_handle(graph->get_handle(snarl->start().node_id(), snarl->start().backward()),
                                                                  graph->get_handle(snarl->end().node_id(), snarl->end().backward()));
                        if (graph->has_edge(deletion_edge)) {
                            edges_to_delete.insert(deletion_edge);
                        }
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
    size_t edge_count = graph->get_edge_count();
    vector<edge_t> edges;
    edges.reserve(edge_count);
    graph->for_each_edge([&](edge_t edge) {
            edges.push_back(edge);
        });
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
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path" << kv.first << endl;
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
                                    SnarlManager& snarl_manager, bool include_endpoints, int64_t min_depth, int64_t min_fragment_len, bool verbose) {

    function<void(function<void(handle_t, const Region*)>)> iterate_handles = [&] (function<void(handle_t, const Region*)> visit_handle) {
        
        visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                                                                                        int64_t start_pos, int64_t end_pos,
                                                                                        bool steps_reversed, const Region* containing_region) {
                                   
                                   pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
                                   for (id_t node_id : contents.first) {
                                       visit_handle(graph->get_handle(node_id), containing_region);
                                   }                                   
                               });
    };

    // todo: duplicating this is very wasteful, and is only happening because edge support was added after the fact.
    // something needs to be refactored in order to fix this, but it's a fairly esoteric codepath and may not be worth it
    function<void(function<void(edge_t, const Region*)>)> iterate_edges = [&] (function<void(edge_t, const Region*)> visit_edge) {
        
        visit_contained_snarls(pp_graph, regions, snarl_manager, include_endpoints, [&](const Snarl* snarl, step_handle_t start_step, step_handle_t end_step,
                                                                                        int64_t start_pos, int64_t end_pos,
                                                                                        bool steps_reversed, const Region* containing_region) {
                                   
                                   pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(snarl, *pp_graph, false);
                                   for (const edge_t& edge : contents.second) {
                                       visit_edge(edge, containing_region);
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
        for (int64_t i = 0; i < context_steps; ++i) {
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

// walk from given edge back to reference path, returning the maximum distance found
static int64_t get_max_deletion(const PathHandleGraph* graph, 
                                const unordered_map<handle_t, vector<int64_t>>& handle_to_position,
                                int64_t context_steps,
                                edge_t edge,
                                const unordered_set<edge_t>& edge_blacklist) {
    
    function<unordered_set<handle_t>(handle_t)> get_ref_context = [&] (handle_t handle) {
        unordered_set<handle_t> context;
        unordered_set<handle_t> ref_context;
        vector<handle_t> cur_handles = {handle};
        for (int64_t i = 0; i < context_steps; ++i) {
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

    // search away from the start of the edge, finding the leftmost and rightmost reference
    // path positions
    unordered_set<handle_t> left_context = get_ref_context(graph->flip(edge.first));
    if (left_context.empty()) {
        return -1;
    }
    int64_t min_pos_left = numeric_limits<int64_t>::max();
    int64_t max_pos_left = -1;
    for (handle_t handle : left_context) {
        auto it = handle_to_position.find(handle);
        assert(it != handle_to_position.end());
        const vector<int64_t>& positions = it->second;
        for (int64_t pos : positions) {
            min_pos_left = min(min_pos_left, pos);
            max_pos_left = max(max_pos_left, pos);
        }
    }

    // search away from the end of the edge, finding the leftmost and rightmost reference
    // path positions
    unordered_set<handle_t> right_context = get_ref_context(edge.second);
    if (right_context.empty()) {
        return -1;
    }
    int64_t min_pos_right = numeric_limits<int64_t>::max();
    int64_t max_pos_right = -1;
    for (handle_t handle : right_context) {
        auto it = handle_to_position.find(handle);
        assert(it != handle_to_position.end());
        const vector<int64_t>& positions = it->second;
        for (int64_t pos : positions) {
            min_pos_right = min(min_pos_right, pos);
            max_pos_right = max(max_pos_right, pos);
        }
    }

    // compute the maximum deletion
    int64_t delta = max(abs(max_pos_left - min_pos_right), abs(max_pos_right - min_pos_left));
    
    return delta;
}

void clip_deletion_edges(MutablePathMutableHandleGraph* graph, int64_t max_deletion,
                         int64_t context_steps,
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

        if (verbose) {
            cerr << "[vg clip]: Found " << length_to_edge.size() << " candidate deletion edges for " << graph->get_path_name(ref_path);
            if (!length_to_edge.empty()) {
                cerr << " with sizes ranging from " << length_to_edge.begin()->first << " to " << length_to_edge.rbegin()->first;
            }
            cerr << endl;
        }

        // for every deletion candidate, walk *back* to the reference path and add it if it forms a deletion
        // we do this one-at-a-time to make sure that the edges are still deletions after the preceding edges were deleted
        int64_t candidate_i = 0;
        for (auto it = length_to_edge.rbegin(); it != length_to_edge.rend(); ++it) {
            const edge_t& edge = it->second;
            assert(graph->has_edge(edge));                                    
            int64_t deletion_size = get_max_deletion(graph, handle_to_position, context_steps, edge, edge_blacklist);

            if (deletion_size > max_deletion) {
                deletion_edges.insert(edge);
                edge_blacklist.insert(edge);
                if (verbose) {
                    if (verbose) {
                        cerr << "[vg clip]: Found deletion edge for candidate " << candidate_i << " on " << graph->get_path_name(ref_path) << ": "
                             << (graph->get_is_reverse(edge.first) ? "<" : ">") << graph->get_id(edge.first)
                             << (graph->get_is_reverse(edge.second) ? "<" : ">") << graph->get_id(edge.second)
                             << " with reference length " << deletion_size << endl;
                    }
                }
            }
            ++candidate_i;
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
            cerr << "[vg-clip]: Creating " << kv.second << " fragments from path" << kv.first << endl;
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


}
