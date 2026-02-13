#include "augref.hpp"
#include <sstream>
#include <algorithm>
#include <queue>
#include <iomanip>

//#define debug

namespace vg {

using namespace std;

const string AugRefCover::augref_suffix = "_alt";

string AugRefCover::make_augref_name(const string& base_path_name, int64_t augref_index) {
    // New naming convention: {base}_{N}_alt
    return base_path_name + "_" + to_string(augref_index) + augref_suffix;
}

bool AugRefCover::is_augref_name(const string& path_name) {
    // Check for pattern: _{digits}_alt at the end
    // Must end with "_alt"
    if (path_name.length() < 6) {  // minimum: "x_1_alt"
        return false;
    }
    if (path_name.substr(path_name.length() - 4) != augref_suffix) {
        return false;
    }
    // Find the underscore before the digits
    size_t alt_pos = path_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = path_name.rfind('_', alt_pos - 1);
    if (underscore_pos == string::npos || underscore_pos == alt_pos - 1) {
        return false;  // no underscore found, or nothing between underscore and _alt
    }
    // Check that everything between underscore and _alt is digits
    for (size_t i = underscore_pos + 1; i < alt_pos; ++i) {
        if (!isdigit(path_name[i])) {
            return false;
        }
    }
    return true;
}

string AugRefCover::parse_base_path(const string& augref_name) {
    if (!is_augref_name(augref_name)) {
        return augref_name;
    }
    // Find _{N}_alt and strip it
    size_t alt_pos = augref_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = augref_name.rfind('_', alt_pos - 1);
    return augref_name.substr(0, underscore_pos);
}

int64_t AugRefCover::parse_augref_index(const string& augref_name) {
    if (!is_augref_name(augref_name)) {
        return -1;
    }
    // Extract N from _{N}_alt
    size_t alt_pos = augref_name.length() - 4;  // position of "_alt"
    size_t underscore_pos = augref_name.rfind('_', alt_pos - 1);
    return stoll(augref_name.substr(underscore_pos + 1, alt_pos - underscore_pos - 1));
}

void AugRefCover::set_augref_sample(const string& sample_name) {
    this->augref_sample_name = sample_name;
}

const string& AugRefCover::get_augref_sample() const {
    return this->augref_sample_name;
}

void AugRefCover::set_verbose(bool verbose) {
    this->verbose = verbose;
}

bool AugRefCover::get_verbose() const {
    return this->verbose;
}

void AugRefCover::clear(MutablePathMutableHandleGraph* graph) {
    vector<path_handle_t> augref_paths_to_remove;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (is_augref_name(graph->get_path_name(path_handle))) {
            augref_paths_to_remove.push_back(path_handle);
        }
    });
    for (path_handle_t path_handle : augref_paths_to_remove) {
        graph->destroy_path(path_handle);
    }
}

void AugRefCover::compute(const PathHandleGraph* graph,
                            SnarlManager* snarl_manager,
                            const unordered_set<path_handle_t>& reference_paths,
                            int64_t minimum_length) {

    // start from scratch
    this->augref_intervals.clear();
    this->interval_snarl_bounds.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        this->augref_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                    graph->path_end(ref_path_handle)));
        this->interval_snarl_bounds.push_back({0, 0});
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (node_to_interval.count(node_id)) {
                cerr << "[augref error]: node " << node_id << " covered by two reference paths,"
                     << " including " << graph->get_path_name(ref_path_handle) << " and "
                     << graph->get_path_name(graph->get_path_handle_of_step(augref_intervals.at(node_to_interval.at(node_id)).first))
                     << ". Augmented reference path support currently requires disjoint acyclic reference paths" << endl;
                exit(1);
            }
            node_to_interval[node_id] = augref_intervals.size() - 1;
        });
    }
    this->num_ref_intervals = this->augref_intervals.size();

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[augref] Selected " << augref_intervals.size() << " rank=0 reference paths" << endl;
#endif

    // we use the path traversal finder for everything
    PathTraversalFinder path_trav_finder(*graph);

    // we collect the augref cover in parallel as a list of path fragments
    size_t thread_count = get_thread_count();
    vector<vector<pair<step_handle_t, step_handle_t>>> augref_intervals_vector(thread_count);
    vector<unordered_map<nid_t, int64_t>> node_to_interval_vector(thread_count);
    vector<vector<pair<nid_t, nid_t>>> snarl_bounds_vector(thread_count);

    // we process top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
        // per-thread output
        vector<pair<step_handle_t, step_handle_t>>& thread_augref_intervals = augref_intervals_vector[omp_get_thread_num()];
        unordered_map<nid_t, int64_t>& thread_node_to_interval = node_to_interval_vector[omp_get_thread_num()];
        vector<pair<nid_t, nid_t>>& thread_snarl_bounds = snarl_bounds_vector[omp_get_thread_num()];

        // capture the top-level snarl boundary node IDs
        nid_t top_snarl_start = snarl->start().node_id();
        nid_t top_snarl_end = snarl->end().node_id();

        vector<const Snarl*> queue = {snarl};

        while(!queue.empty()) {
            const Snarl* cur_snarl = queue.back();
            queue.pop_back();

            // get the snarl cover
            compute_snarl(*cur_snarl, path_trav_finder, minimum_length,
                          thread_augref_intervals,
                          thread_node_to_interval,
                          top_snarl_start, top_snarl_end,
                          thread_snarl_bounds);

            // recurse on the children
            const vector<const Snarl*>& children = snarl_manager->children_of(cur_snarl);
            for (const Snarl* child_snarl : children) {
                queue.push_back(child_snarl);
            }
        }
    });

    // now we need to fold up the thread covers
    for (int64_t t = 0; t < thread_count; ++t) {
#ifdef debug
#pragma omp critical(cerr)
        cerr << "Adding " << augref_intervals_vector[t].size() << " augref intervals from thread " << t << endl;
#endif
        // important to go through function rather than do a raw copy since
        // inter-top-level snarl merging may need to happen
        for (int64_t j = 0; j < augref_intervals_vector[t].size(); ++j) {
            // the true flag at the end disables the overlap check. since they were computed
            // in separate threads, snarls can overlap by a single node
            const pair<step_handle_t, step_handle_t>& interval = augref_intervals_vector[t][j];
            if (interval.first != graph->path_end(graph->get_path_handle_of_step(interval.first))) {
                add_interval(this->augref_intervals, this->node_to_interval, augref_intervals_vector[t][j], true,
                             &this->interval_snarl_bounds, snarl_bounds_vector[t][j]);
            }
        }
        augref_intervals_vector[t].clear();
        node_to_interval_vector[t].clear();
        snarl_bounds_vector[t].clear();
    }

    // second pass: greedily cover any nodes not covered by snarl traversals
    fill_uncovered_nodes(minimum_length);

    // debug: verify all nodes are covered
    verify_cover();

    // remove any intervals that were made redundant by add_interval
    defragment_intervals();
}

void AugRefCover::fill_uncovered_nodes(int64_t minimum_length) {
    // Collect all uncovered nodes and the paths that pass through them
    unordered_set<nid_t> uncovered_nodes;
    map<string, path_handle_t> candidate_paths;  // sorted by name for deterministic ordering

    graph->for_each_handle([&](handle_t handle) {
        nid_t node_id = graph->get_id(handle);
        if (!node_to_interval.count(node_id)) {
            // Node is not covered - find paths through it
            uncovered_nodes.insert(node_id);
            graph->for_each_step_on_handle(handle, [&](step_handle_t step) {
                path_handle_t path_handle = graph->get_path_handle_of_step(step);
                string path_name = graph->get_path_name(path_handle);
                // Skip existing augref paths
                if (!is_augref_name(path_name)) {
                    candidate_paths[path_name] = path_handle;
                }
                return true;
            });
        }
    });

    if (uncovered_nodes.empty()) {
        return;
    }

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[augref] fill_uncovered_nodes: " << uncovered_nodes.size() << " uncovered nodes, "
         << candidate_paths.size() << " candidate paths" << endl;
#endif

    // Greedily walk through paths, creating intervals for contiguous uncovered sequences
    for (const auto& name_path : candidate_paths) {
        if (uncovered_nodes.empty()) {
            break;
        }

        bool in_interval = false;
        step_handle_t interval_start;
        step_handle_t interval_end;
        int64_t interval_length = 0;
        unordered_set<nid_t> interval_nodes;  // track nodes in current interval for cycle detection

        // Helper to close the current interval and add it if long enough
        auto close_interval = [&]() {
            if (in_interval) {
                if (interval_length >= minimum_length) {
                    add_interval(this->augref_intervals, this->node_to_interval,
                                 make_pair(interval_start, interval_end), true,
                                 &this->interval_snarl_bounds, {0, 0});
                    for (nid_t nid : interval_nodes) {
                        uncovered_nodes.erase(nid);
                    }
                }
                in_interval = false;
            }
        };

        graph->for_each_step_in_path(name_path.second, [&](step_handle_t step) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));

            if (uncovered_nodes.count(node_id)) {
                if (interval_nodes.count(node_id)) {
                    // This node is already in the current interval — close to avoid cycle,
                    // then start a new interval at this node (same logic as get_uncovered_intervals)
                    close_interval();
                }
                if (!in_interval) {
                    // Start a new interval
                    in_interval = true;
                    interval_start = step;
                    interval_length = 0;
                    interval_nodes.clear();
                }
                interval_end = graph->get_next_step(step);
                interval_length += graph->get_length(graph->get_handle_of_step(step));
                interval_nodes.insert(node_id);
            } else {
                // This node is already covered — close current interval
                close_interval();
            }
            return true;
        });

        // Don't forget to close any interval at the end of the path
        close_interval();
    }

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[augref] fill_uncovered_nodes: " << uncovered_nodes.size() << " nodes still uncovered after second pass" << endl;
#endif
}

void AugRefCover::load(const PathHandleGraph* graph,
                         const unordered_set<path_handle_t>& reference_paths) {
    // start from scratch
    this->augref_intervals.clear();
    this->interval_snarl_bounds.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (graph->get_is_reverse(graph->get_handle_of_step(step_handle))) {
                cerr << "[augref] error: Reversed step " << node_id << " found in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All augref fragments must be forward-only." << endl;
                exit(1);
            }
            if (node_to_interval.count(node_id)) {
                cerr << "[augref] error: Cycle found on node " << node_id << " in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All augref fragments must be acyclic." << endl;
                exit(1);
            }
            node_to_interval[node_id] = augref_intervals.size();
        });
        this->augref_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                    graph->path_end(ref_path_handle)));
        this->interval_snarl_bounds.push_back({0, 0});
    }
    this->num_ref_intervals = this->augref_intervals.size();

    // load existing augref paths from the graph
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (is_augref_name(path_name)) {
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                node_to_interval[graph->get_id(graph->get_handle_of_step(step_handle))] = augref_intervals.size();
            });
            this->augref_intervals.push_back(make_pair(graph->path_begin(path_handle),
                                                        graph->path_end(path_handle)));
            this->interval_snarl_bounds.push_back({0, 0});
        }
    });
}

void AugRefCover::apply(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));
#ifdef debug
    cerr << "applying augref cover with " << this->num_ref_intervals << " ref intervals "
         << " and " << this->augref_intervals.size() << " total intervals" << endl;
#endif

    // If augref_sample_name is set, first copy base reference paths to the new sample
    if (!augref_sample_name.empty()) {
        // Collect reference path handles from the reference intervals
        unordered_set<path_handle_t> reference_paths;
        for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
            reference_paths.insert(graph->get_path_handle_of_step(augref_intervals[i].first));
        }
        copy_base_paths_to_sample(mutable_graph, reference_paths);
    }

    // Reset augref counters for each base path
    base_path_augref_counter.clear();

    // First pass: determine the maximum existing augref index for each base path
    // This ensures we don't overwrite existing augref paths
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_augref_name(path_name)) {
            string base = parse_base_path(path_name);
            int64_t idx = parse_augref_index(path_name);
            if (base_path_augref_counter.count(base)) {
                base_path_augref_counter[base] = max(base_path_augref_counter[base], idx);
            } else {
                base_path_augref_counter[base] = idx;
            }
        }
    });

    // write the augref paths
    int64_t written_intervals = 0;
    int64_t written_length = 0;
    int64_t skipped_intervals = 0;
    for (int64_t i = this->num_ref_intervals; i < this->augref_intervals.size(); ++i) {
        // Skip empty intervals (these can be created by defragment_intervals or merging)
        path_handle_t interval_path = graph->get_path_handle_of_step(augref_intervals[i].first);
        if (augref_intervals[i].first == graph->path_end(interval_path)) {
            skipped_intervals++;
            continue;
        }

        // Find the reference path this augref path extends from by tracing back to reference
        nid_t first_node = graph->get_id(graph->get_handle_of_step(augref_intervals[i].first));
        vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, true);

        // Get the reference path name from the reference node
        string base_path_name;
        if (!ref_nodes.empty()) {
            nid_t ref_node_id = ref_nodes.at(0).second;
            int64_t ref_interval_idx = this->node_to_interval.at(ref_node_id);
            path_handle_t ref_path_handle = graph->get_path_handle_of_step(augref_intervals[ref_interval_idx].first);
            base_path_name = graph->get_path_name(ref_path_handle);
            // Strip any subrange from the reference path name
            subrange_t subrange;
            base_path_name = Paths::strip_subrange(base_path_name, &subrange);
        } else {
            // Fallback to source path if no reference found (shouldn't happen)
            path_handle_t source_path_handle = mutable_graph->get_path_handle_of_step(augref_intervals[i].first);
            base_path_name = graph->get_path_name(source_path_handle);
            subrange_t subrange;
            base_path_name = Paths::strip_subrange(base_path_name, &subrange);
        }

        // If augref_sample_name is set, replace the sample in base_path_name
        if (!augref_sample_name.empty()) {
            PathSense sense;
            string sample, locus;
            size_t haplotype, phase_block;
            subrange_t subrange;
            PathMetadata::parse_path_name(base_path_name, sense, sample, locus, haplotype, phase_block, subrange);

            if (sample.empty()) {
                // Simple path name - prepend augref sample
                base_path_name = augref_sample_name + "#0#" + base_path_name;
            } else {
                // Replace sample with augref sample
                base_path_name = PathMetadata::create_path_name(sense, augref_sample_name, locus, haplotype, phase_block, subrange);
            }
        }

        // Get next available augref index for this base path
        int64_t augref_index = ++base_path_augref_counter[base_path_name];

        // Create the augref path name
        string augref_name = make_augref_name(base_path_name, augref_index);

        // Create the path as REFERENCE sense
        path_handle_t augref_handle = mutable_graph->create_path_handle(augref_name, false);

        int64_t interval_length = 0;
        for (step_handle_t step_handle = augref_intervals[i].first; step_handle != augref_intervals[i].second;
             step_handle = mutable_graph->get_next_step(step_handle)) {
            mutable_graph->append_step(augref_handle, mutable_graph->get_handle_of_step(step_handle));
            interval_length += mutable_graph->get_length(mutable_graph->get_handle_of_step(step_handle));
        }
        written_intervals++;
        written_length += interval_length;
    }

#ifdef debug
    cerr << "[augref] apply: wrote " << written_intervals << " augref paths (" << written_length << " bp), skipped " << skipped_intervals << " empty intervals" << endl;
#endif

    this->forwardize_augref_paths(mutable_graph);
}

int64_t AugRefCover::get_rank(nid_t node_id) const {
    // search back to reference in order to find the rank.
    vector<pair<int64_t, nid_t>> ref_steps = this->get_reference_nodes(node_id, true);
    // Return -1 if node is in a disconnected component that can't reach reference
    return ref_steps.empty() ? -1 : ref_steps.at(0).first;
}

const vector<pair<step_handle_t, step_handle_t>>& AugRefCover::get_intervals() const {
    return this->augref_intervals;
}

const pair<step_handle_t, step_handle_t>* AugRefCover::get_interval(nid_t node_id) const {
    if (this->node_to_interval.count(node_id)) {
        return &this->augref_intervals.at(node_to_interval.at(node_id));
    }
    return nullptr;
}

int64_t AugRefCover::get_num_ref_intervals() const {
    return this->num_ref_intervals;
}

void AugRefCover::compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                                  vector<pair<step_handle_t, step_handle_t>>& thread_augref_intervals,
                                  unordered_map<nid_t, int64_t>& thread_node_to_interval,
                                  nid_t top_snarl_start, nid_t top_snarl_end,
                                  vector<pair<nid_t, nid_t>>& thread_snarl_bounds) {

    // start by finding the path traversals through the snarl
    vector<vector<step_handle_t>> travs;
    vector<string> trav_names;
    {
        pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs = path_trav_finder.find_path_traversals(snarl);
        travs.reserve(path_travs.first.size());

        // reduce protobuf usage by going back to vector of steps instead of keeping SnarlTraversals around
        for (int64_t i = 0; i < path_travs.first.size(); ++i) {
            string trav_path_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
            if (is_augref_name(trav_path_name)) {
                // we ignore existing (off-reference) augref paths
#ifdef debug
                cerr << "Warning : skipping existing augref traversal " << trav_path_name << endl;
#endif
                continue;
            }
            bool reversed = false;
            if (graph->get_is_reverse(graph->get_handle_of_step(path_travs.second[i].first)) != snarl.start().backward()) {
                reversed = true;
            }
            assert((graph->get_is_reverse(graph->get_handle_of_step(path_travs.second[i].second)) != snarl.end().backward()) == reversed);
            vector<step_handle_t> trav;
            trav.reserve(path_travs.first[i].visit_size());
            bool done = false;
            function<step_handle_t(step_handle_t)> visit_next_step = [&](step_handle_t step_handle) {
                return reversed ? graph->get_previous_step(step_handle) : graph->get_next_step(step_handle);
            };
            for (step_handle_t step_handle = path_travs.second[i].first; !done; step_handle = visit_next_step(step_handle)) {
                trav.push_back(step_handle);
                if (step_handle == path_travs.second[i].second) {
                    done = true;
                }
            }
            if (reversed) {
                std::reverse(trav.begin(), trav.end());
            }
            travs.push_back(trav);
            trav_names.push_back(trav_path_name);
        }
    }
#ifdef debug
#pragma omp critical(cerr)
    cerr << "doing snarl " << pb2json(snarl.start()) << "-" << pb2json(snarl.end()) << " with " << travs.size() << " travs" << endl;
#endif

    // build an initial ranked list of candidate traversal fragments
    vector<RankedFragment> ranked_trav_fragments;
    for (int64_t trav_idx = 0; trav_idx < travs.size(); ++trav_idx) {
        // only a reference traversal (or deletion that we don't need to consider)
        // will have its first two nodes covered
        if (this->node_to_interval.count(graph->get_id(graph->get_handle_of_step(travs[trav_idx][0]))) &&
            this->node_to_interval.count(graph->get_id(graph->get_handle_of_step(travs[trav_idx][1])))) {
            continue;
        }

        const vector<step_handle_t>& trav = travs.at(trav_idx);
        vector<pair<int64_t, int64_t>> uncovered_intervals = get_uncovered_intervals(trav, thread_node_to_interval);

        for (const auto& uncovered_interval : uncovered_intervals) {
            unordered_set<nid_t> cycle_check;
            bool cyclic = false;
            int64_t interval_length = 0;
            for (int64_t i = uncovered_interval.first; i < uncovered_interval.second && !cyclic; ++i) {
                handle_t handle = graph->get_handle_of_step(trav[i]);
                interval_length += graph->get_length(handle);
                nid_t node_id = graph->get_id(handle);
                if (cycle_check.count(node_id)) {
                    cyclic = true;
                } else {
                    cycle_check.insert(node_id);
                }
            }
            if (!cyclic && interval_length >= minimum_length) {
                int64_t trav_coverage = get_coverage(trav, uncovered_interval);
                ranked_trav_fragments.push_back({trav_coverage, &trav_names[trav_idx], trav_idx, uncovered_interval});
            }
        }
    }

    // put the fragments into a max heap
    std::make_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end());

    // now greedily pull out traversal intervals from the ranked list until none are left
    while (!ranked_trav_fragments.empty()) {

        // get the best scoring (max) fragment from heap
        auto best_stats_fragment = ranked_trav_fragments.front();
        std::pop_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end());
        ranked_trav_fragments.pop_back();

        const vector<step_handle_t>& trav = travs.at(best_stats_fragment.trav_idx);
        const pair<int64_t, int64_t>& uncovered_interval = best_stats_fragment.fragment;

#ifdef debug
#pragma omp critical(cerr)
        {
        cerr << "best trav: ";
        for (auto xx : trav) cerr << " " << graph->get_id(graph->get_handle_of_step(xx));
        cerr << endl << "uncovered interval [" << uncovered_interval.first << "," << uncovered_interval.second << "]" <<endl;
        }
#endif

        // our traversal may have been partially covered by a different iteration, if so, we need to break it up
        // and continue
        vector<pair<int64_t, int64_t>> chopped_intervals;
        int64_t cur_start = -1;
        bool chopped = false;
        for (int64_t i = uncovered_interval.first; i < uncovered_interval.second; ++i) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(trav[i]));
            bool covered = this->node_to_interval.count(node_id) || thread_node_to_interval.count(node_id);
            if (!covered && cur_start == -1) {
                cur_start = i;
            } else if (covered) {
                chopped = true;
                if (cur_start != -1) {
                    chopped_intervals.push_back(make_pair(cur_start, i));
                    cur_start = -1;
                }
            }
        }
        if (cur_start != -1) {
            chopped_intervals.push_back(make_pair(cur_start, uncovered_interval.second));
        }
        if (chopped) {
            for (const pair<int64_t, int64_t>& chopped_interval : chopped_intervals) {
                int64_t chopped_trav_length = 0;
                for (int64_t i = chopped_interval.first; i < chopped_interval.second; ++i) {
                    chopped_trav_length += graph->get_length(graph->get_handle_of_step(trav[i]));
                }
                if (chopped_trav_length >= minimum_length) {
                    int64_t trav_coverage = get_coverage(trav, chopped_interval);
                    ranked_trav_fragments.push_back({trav_coverage, best_stats_fragment.name, best_stats_fragment.trav_idx, chopped_interval});
                    std::push_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end());
                }
            }
            continue;
        }
        pair<step_handle_t, step_handle_t> new_interval = make_pair(trav.at(uncovered_interval.first),
                                                                    graph->get_next_step(trav.at(uncovered_interval.second - 1)));
#ifdef debug
        int64_t interval_length = uncovered_interval.second - uncovered_interval.first;
#pragma omp critical(cerr)
        cerr << "adding interval with length " << interval_length << endl;
#endif
        add_interval(thread_augref_intervals, thread_node_to_interval, new_interval, false,
                     &thread_snarl_bounds, {top_snarl_start, top_snarl_end});
    }
}

vector<pair<int64_t, int64_t>> AugRefCover::get_uncovered_intervals(const vector<step_handle_t>& trav,
                                                                      const unordered_map<nid_t, int64_t>& thread_node_to_interval) {

    vector<pair<int64_t, int64_t>> intervals;
    int64_t start = -1;
    unordered_set<nid_t> dupe_check;
    for (size_t i = 0; i < trav.size(); ++i) {
        nid_t node_id = graph->get_id(graph->get_handle_of_step(trav[i]));
        bool covered = this->node_to_interval.count(node_id) || thread_node_to_interval.count(node_id);
        // we break at dupes even if uncovered -- never want same id twice in an interval
        bool dupe = !covered && dupe_check.count(node_id);
        dupe_check.insert(node_id);
        if (covered || dupe) {
            if (start != -1) {
                intervals.push_back(make_pair(start, i));
            }
            start = dupe ? i : -1;
        } else {
            if (start == -1) {
                start = i;
            }
        }
    }
    if (start != -1) {
        intervals.push_back(make_pair(start, trav.size()));
    }
    return intervals;
}

bool AugRefCover::add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_augref_intervals,
                                 unordered_map<nid_t, int64_t>& thread_node_to_interval,
                                 const pair<step_handle_t, step_handle_t>& new_interval,
                                 bool global,
                                 vector<pair<nid_t, nid_t>>* snarl_bounds_vec,
                                 pair<nid_t, nid_t> snarl_bounds) {

#ifdef debug
#pragma omp critical(cerr)
    cerr << "adding interval " << graph->get_path_name(graph->get_path_handle_of_step(new_interval.first))
         << (graph->get_is_reverse(graph->get_handle_of_step(new_interval.first)) ? "<" : ">")
         << ":" << graph->get_id(graph->get_handle_of_step(new_interval.first));
    if (new_interval.second == graph->path_end(graph->get_path_handle_of_step(new_interval.first))) {
        cerr << "PATH_END" << endl;
    }  else {
        cerr << "-" << (graph->get_is_reverse(graph->get_handle_of_step(new_interval.second)) ? "<" : ">")
             << graph->get_id(graph->get_handle_of_step(new_interval.second)) << endl;
    }
#endif
    bool merged = false;
    int64_t merged_interval_idx = -1;
    path_handle_t path_handle = graph->get_path_handle_of_step(new_interval.first);

    // check the before-first step. if it's in an interval then it must be immediately
    // preceeding so we merge the new interval to the end of the found interval
    step_handle_t before_first_step = graph->get_previous_step(new_interval.first);
    if (before_first_step != graph->path_front_end(graph->get_path_handle_of_step(before_first_step))) {
        nid_t prev_node_id = graph->get_id(graph->get_handle_of_step(before_first_step));
        if (thread_node_to_interval.count(prev_node_id)) {
            int64_t prev_idx = thread_node_to_interval[prev_node_id];
            pair<step_handle_t, step_handle_t>& prev_interval = thread_augref_intervals[prev_idx];
            if (graph->get_path_handle_of_step(prev_interval.first) == path_handle) {
                if (prev_interval.second == new_interval.first ||
                    (global && graph->get_previous_step(prev_interval.second) == new_interval.first)) {
#ifdef debug
#pragma omp critical(cerr)
                    cerr << "prev interval found" << graph->get_path_name(graph->get_path_handle_of_step(prev_interval.first))
                         << ":" << (graph->get_is_reverse(graph->get_handle_of_step(prev_interval.first)) ? "<" : ">")
                         << graph->get_id(graph->get_handle_of_step(prev_interval.first));
                    if (prev_interval.second == graph->path_end(graph->get_path_handle_of_step(prev_interval.first))) {
                        cerr << "PATH_END" << endl;
                    } else {
                         cerr << "-" << (graph->get_is_reverse(graph->get_handle_of_step(prev_interval.second)) ? "<" : ">")
                              << graph->get_id(graph->get_handle_of_step(prev_interval.second)) << endl;
                    }
#endif
                    prev_interval.second = new_interval.second;
                    merged = true;
                    merged_interval_idx = prev_idx;
                }
            }
        }
    }

    // check the end step. if it's in an interval then it must be immediately
    // following we merge the new interval to the front of the found interval
    int64_t deleted_idx = -1;
    step_handle_t deleted_interval_first, deleted_interval_second;  // save before decommission
    if (new_interval.second != graph->path_end(graph->get_path_handle_of_step(new_interval.second))) {
        nid_t next_node_id = graph->get_id(graph->get_handle_of_step(new_interval.second));
        if (thread_node_to_interval.count(next_node_id)) {
            int64_t next_idx = thread_node_to_interval[next_node_id];
            pair<step_handle_t, step_handle_t>& next_interval = thread_augref_intervals[next_idx];
            path_handle_t next_path = graph->get_path_handle_of_step(next_interval.first);
            if (graph->get_path_handle_of_step(next_interval.first) == path_handle) {
                if (next_interval.first == new_interval.second ||
                    (global && next_interval.first == graph->get_previous_step(new_interval.second))) {
#ifdef debug
#pragma omp critical(cerr)
                    cerr << "next interval found" << graph->get_path_name(graph->get_path_handle_of_step(next_interval.first))
                         << ":" << graph->get_id(graph->get_handle_of_step(next_interval.first));
                    if (next_interval.second == graph->path_end(graph->get_path_handle_of_step(next_interval.second))) {
                        cerr << "PATH_END" << endl;
                    } else {
                         cerr << "-" << graph->get_id(graph->get_handle_of_step(next_interval.second)) << endl;
                    }
#endif
                    if (merged == true) {
                        // save the interval bounds BEFORE decommissioning
                        deleted_interval_first = next_interval.first;
                        deleted_interval_second = next_interval.second;
                        deleted_idx = next_idx;
                        // decomission next_interval
                        next_interval.first = graph->path_end(next_path);
                        next_interval.second = graph->path_front_end(next_path);
                        if (snarl_bounds_vec) {
                            (*snarl_bounds_vec)[next_idx] = {0, 0};
                        }
                        // extend the previous interval right to cover both new_interval and the deleted next_interval
                        thread_augref_intervals[merged_interval_idx].second = deleted_interval_second;
                    } else {
                        // extend next_interval left
                        next_interval.first = new_interval.first;
                        merged = true;
                        merged_interval_idx = next_idx;
                    }
                }
            }
        }
    }

    // add the interval to the local (thread safe) structures
    if (!merged) {
        merged_interval_idx = thread_augref_intervals.size();
        thread_augref_intervals.push_back(new_interval);
        if (snarl_bounds_vec) {
            snarl_bounds_vec->push_back(snarl_bounds);
        }
    }
    for (step_handle_t step = new_interval.first; step != new_interval.second; step = graph->get_next_step(step)) {
        thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = merged_interval_idx;
    }
    if (deleted_idx >= 0) {
        // move the links to the deleted interval to the merged interval
        // use saved bounds since the interval has been decommissioned
        for (step_handle_t step = deleted_interval_first; step != deleted_interval_second; step = graph->get_next_step(step)) {
            thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = merged_interval_idx;
        }
    }
    return !merged;
}

void AugRefCover::defragment_intervals() {
    vector<pair<step_handle_t, step_handle_t>> new_intervals;
    vector<pair<nid_t, nid_t>> new_snarl_bounds;
    this->node_to_interval.clear();
    for (int64_t i = 0; i < this->augref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->augref_intervals[i];
        path_handle_t path_handle = graph->get_path_handle_of_step(interval.first);
        if (interval.first != graph->path_end(path_handle)) {
            new_intervals.push_back(interval);
            new_snarl_bounds.push_back(this->interval_snarl_bounds[i]);
        }
    }
    this->augref_intervals = std::move(new_intervals);
    this->interval_snarl_bounds = std::move(new_snarl_bounds);
    for (int64_t i = 0; i < this->augref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->augref_intervals[i];
        for (step_handle_t step = interval.first; step != interval.second; step = graph->get_next_step(step)) {
            this->node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = i;
        }
    }
}

int64_t AugRefCover::get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval) {
    int64_t coverage = 0;

    for (int64_t i = uncovered_interval.first; i < uncovered_interval.second; ++i) {
        const step_handle_t& step = trav[i];
        handle_t handle = graph->get_handle_of_step(step);
        vector<step_handle_t> all_steps = graph->steps_of_handle(handle);
        int64_t length = graph->get_length(handle);
        coverage += length * all_steps.size();
    }

    return coverage;
}

// Ensure all nodes in augref paths are in forward orientation
void AugRefCover::forwardize_augref_paths(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));

    unordered_map<nid_t, nid_t> id_map;
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_augref_name(path_name)) {
            size_t fw_count = 0;
            size_t total_steps = 0;
            mutable_graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = mutable_graph->get_handle_of_step(step_handle);
                if (mutable_graph->get_is_reverse(handle)) {
                    handle_t flipped_handle = mutable_graph->create_handle(mutable_graph->get_sequence(handle));
                    id_map[mutable_graph->get_id(flipped_handle)] = mutable_graph->get_id(handle);
                    mutable_graph->follow_edges(handle, true, [&](handle_t prev_handle) {
                        if (mutable_graph->get_id(prev_handle) != mutable_graph->get_id(handle)) {
                            mutable_graph->create_edge(prev_handle, flipped_handle);
                        }
                    });
                    mutable_graph->follow_edges(handle, false, [&](handle_t next_handle) {
                        if (mutable_graph->get_id(handle) != mutable_graph->get_id(next_handle)) {
                            mutable_graph->create_edge(flipped_handle, next_handle);
                        }
                    });
                    // self-loop cases we punted on above:
                    if (mutable_graph->has_edge(handle, handle)) {
                        mutable_graph->create_edge(flipped_handle, flipped_handle);
                    }
                    if (mutable_graph->has_edge(handle, mutable_graph->flip(handle))) {
                        mutable_graph->create_edge(flipped_handle, mutable_graph->flip(flipped_handle));
                    }
                    if (mutable_graph->has_edge(mutable_graph->flip(handle), handle)) {
                        mutable_graph->create_edge(mutable_graph->flip(flipped_handle), flipped_handle);
                    }
                    vector<step_handle_t> steps = mutable_graph->steps_of_handle(handle);
                    size_t ref_count = 0;
                    for (step_handle_t step : steps) {
                        if (mutable_graph->get_path_handle_of_step(step) == path_handle) {
                            ++ref_count;
                        }
                        step_handle_t next_step = mutable_graph->get_next_step(step);
                        handle_t new_handle = mutable_graph->get_is_reverse(mutable_graph->get_handle_of_step(step)) ? flipped_handle :
                            mutable_graph->flip(flipped_handle);
                        mutable_graph->rewrite_segment(step, next_step, {new_handle});
                    }
                    if (ref_count > 1) {
                        cerr << "[augref] error: Cycle detected in augref path " << path_name << " at node " << mutable_graph->get_id(handle) << endl;
                        exit(1);
                    }
                    ++fw_count;
                    assert(mutable_graph->steps_of_handle(handle).empty());
                    dynamic_cast<DeletableHandleGraph*>(mutable_graph)->destroy_handle(handle);
                }
                ++total_steps;
            });
        }
    });

    // rename all the ids back to what they were (so nodes keep their ids, just get flipped around)
    mutable_graph->reassign_node_ids([&id_map](nid_t new_id) {
        return id_map.count(new_id) ? id_map[new_id] : new_id;
    });

    // do a check just to be sure
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_augref_name(path_name)) {
            mutable_graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = mutable_graph->get_handle_of_step(step_handle);
                if (mutable_graph->get_is_reverse(handle)) {
                    cerr << "[augref] error: Failed to fowardize node " << mutable_graph->get_id(handle) << " in path " << path_name << endl;
                    exit(1);
                }
            });
        }
    });
}

vector<pair<int64_t, nid_t>> AugRefCover::get_reference_nodes(nid_t node_id, bool first) const {

    // search back to reference in order to find the rank.
    unordered_set<nid_t> visited;
    priority_queue<pair<int64_t, nid_t>> queue;
    queue.push(make_pair(0, node_id));

    nid_t current_id;
    int64_t distance = 0;

    // output reference intervals
    vector<pair<int64_t, nid_t>> output_reference_nodes;

    while (!queue.empty()) {
        std::tie(distance, current_id) = queue.top();
        queue.pop();

        if (!visited.count(current_id)) {

            visited.insert(current_id);

            if (this->node_to_interval.count(current_id)) {
                int64_t interval_idx = this->node_to_interval.at(current_id);

                const pair<step_handle_t, step_handle_t>& augref_interval = this->augref_intervals.at(interval_idx);

                // we've hit the reference, fish out its step and stop searching.
                if (interval_idx < this->num_ref_intervals) {
                    output_reference_nodes.push_back(make_pair(distance, current_id));
                    if (first) {
                        break;
                    }
                    continue;
                }

                // search out of the snarl -- any parent traversals will overlap here
                graph->follow_edges(graph->get_handle_of_step(augref_interval.first), true, [&](handle_t prev) {
                    queue.push(make_pair(distance + 1, graph->get_id(prev)));
                });
                // hack around gbwtgraph bug (feature?) that does not let you decrement path_end
                path_handle_t path_handle = graph->get_path_handle_of_step(augref_interval.first);
                step_handle_t last_step;
                if (augref_interval.second == graph->path_end(path_handle)) {
                    last_step = graph->path_back(path_handle);
                } else {
                    last_step = graph->get_previous_step(augref_interval.second);
                }
                graph->follow_edges(graph->get_handle_of_step(last_step), false, [&](handle_t next) {
                    queue.push(make_pair(distance + 1, graph->get_id(next)));
                });

            } else {
                // revert to graph search if node not in interval (distance doesn't increase -- we only count intervals)
                graph->follow_edges(graph->get_handle(current_id), false, [&](handle_t next) {
                    queue.push(make_pair(distance, graph->get_id(next)));
                });
                graph->follow_edges(graph->get_handle(current_id), true, [&](handle_t next) {
                    queue.push(make_pair(distance, graph->get_id(next)));
                });
            }

        }
    }

    // Note: output_reference_nodes may be empty if the node is in a disconnected
    // component that cannot trace back to any reference interval (e.g., after clipping)
    return output_reference_nodes;
}

void AugRefCover::verify_cover() const {
    if (!verbose) {
        return;
    }

    // Check that every node in the graph is covered by the augref cover
    int64_t total_nodes = 0;
    int64_t total_length = 0;
    int64_t ref_nodes = 0;
    int64_t ref_length = 0;
    int64_t alt_nodes = 0;
    int64_t alt_length = 0;
    int64_t uncovered_nodes = 0;
    int64_t uncovered_length = 0;

    graph->for_each_handle([&](handle_t handle) {
        nid_t node_id = graph->get_id(handle);
        int64_t node_len = graph->get_length(handle);
        total_nodes++;
        total_length += node_len;

        if (!node_to_interval.count(node_id)) {
            // Node is not covered (expected when min-augref-len > 0)
            uncovered_nodes++;
            uncovered_length += node_len;
        } else {
            int64_t interval_idx = node_to_interval.at(node_id);
            if (interval_idx < num_ref_intervals) {
                ref_nodes++;
                ref_length += node_len;
            } else {
                alt_nodes++;
                alt_length += node_len;
            }
        }
    });

    cerr << "[augref] verify_cover summary:" << endl
         << "  Total nodes: " << total_nodes << " (" << total_length << " bp)" << endl
         << "  Reference nodes: " << ref_nodes << " (" << ref_length << " bp)" << endl
         << "  Alt nodes: " << alt_nodes << " (" << alt_length << " bp)" << endl
         << "  Uncovered nodes: " << uncovered_nodes << " (" << uncovered_length << " bp)" << endl
         << "  Intervals: " << num_ref_intervals << " ref + " << (augref_intervals.size() - num_ref_intervals) << " alt" << endl;
}

void AugRefCover::copy_base_paths_to_sample(MutablePathMutableHandleGraph* mutable_graph,
                                              const unordered_set<path_handle_t>& reference_paths) {
    if (augref_sample_name.empty()) {
        return;
    }

    for (const path_handle_t& ref_path : reference_paths) {
        string ref_name = mutable_graph->get_path_name(ref_path);

        // Parse the path name to extract components
        // rGFA format: SAMPLE#HAPLOTYPE#CONTIG or just CONTIG
        PathSense sense;
        string sample, locus;
        size_t haplotype, phase_block;
        subrange_t subrange;
        PathMetadata::parse_path_name(ref_name, sense, sample, locus, haplotype, phase_block, subrange);

        // Create new path name with augref sample
        string new_name;
        if (sample.empty()) {
            // Simple path name (no sample info) - prepend augref sample
            new_name = augref_sample_name + "#0#" + ref_name;
        } else {
            // Replace sample with augref sample
            new_name = PathMetadata::create_path_name(sense, augref_sample_name, locus, haplotype, phase_block, subrange);
        }

        // Check if path already exists
        if (mutable_graph->has_path(new_name)) {
#ifdef debug
            cerr << "[augref] copy_base_paths_to_sample: path " << new_name << " already exists, skipping" << endl;
#endif
            continue;
        }

        // Create the new path with same sense as original
        path_handle_t new_path = mutable_graph->create_path_handle(new_name, false);

        // Copy all steps from original path to new path
        mutable_graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
            mutable_graph->append_step(new_path, mutable_graph->get_handle_of_step(step));
            return true;
        });

#ifdef debug
        cerr << "[augref] copy_base_paths_to_sample: copied " << ref_name << " -> " << new_name << endl;
#endif
    }
}

void AugRefCover::write_augref_segments(ostream& os) {
    // Track augref counters to predict path names (same logic as apply())
    unordered_map<string, int64_t> local_augref_counter;

    // First pass: find maximum existing augref index for each base path
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (is_augref_name(path_name)) {
            string base = parse_base_path(path_name);
            int64_t idx = parse_augref_index(path_name);
            if (local_augref_counter.count(base)) {
                local_augref_counter[base] = max(local_augref_counter[base], idx);
            } else {
                local_augref_counter[base] = idx;
            }
        }
    });

    // Pre-compute reference node positions by walking all reference intervals once.
    // Maps node ID -> (ref_path_handle, offset of node start on ref path)
    unordered_map<nid_t, pair<path_handle_t, int64_t>> ref_node_positions;
    for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
        const pair<step_handle_t, step_handle_t>& ref_interval = this->augref_intervals[i];
        path_handle_t ref_path_handle = graph->get_path_handle_of_step(ref_interval.first);
        int64_t offset = 0;
        for (step_handle_t step = ref_interval.first; step != ref_interval.second;
             step = graph->get_next_step(step)) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step));
            ref_node_positions[node_id] = make_pair(ref_path_handle, offset);
            offset += graph->get_length(graph->get_handle_of_step(step));
        }
    }

    // Collect the specific steps whose offsets we need (interval boundaries only).
    // This avoids caching every step of every source path.
    unordered_set<step_handle_t> needed_steps;
    unordered_set<path_handle_t> source_paths_needed;
    for (int64_t i = this->num_ref_intervals; i < this->augref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->augref_intervals[i];
        path_handle_t ph = graph->get_path_handle_of_step(interval.first);
        source_paths_needed.insert(ph);
        needed_steps.insert(interval.first);
        needed_steps.insert(interval.second);
    }
    // Walk each source path once, caching offsets only for needed steps.
    // Also cache path_end -> total path length for computing interval end offsets.
    unordered_map<step_handle_t, int64_t> step_offset_cache;
    for (const path_handle_t& ph : source_paths_needed) {
        int64_t offset = 0;
        for (step_handle_t step = graph->path_begin(ph);
             step != graph->path_end(ph);
             step = graph->get_next_step(step)) {
            if (needed_steps.count(step)) {
                step_offset_cache[step] = offset;
            }
            offset += graph->get_length(graph->get_handle_of_step(step));
        }
        // Cache path_end sentinel with total path length
        step_handle_t path_end = graph->path_end(ph);
        if (needed_steps.count(path_end)) {
            step_offset_cache[path_end] = offset;
        }
    }

    // Cache resolved source display names per path handle to avoid
    // repeated parse_path_name/create_path_name per interval.
    // Stores (display_name, subrange_offset) per source path.
    unordered_map<path_handle_t, pair<string, int64_t>> source_display_cache;

    // Write each augref interval
    for (int64_t i = this->num_ref_intervals; i < this->augref_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->augref_intervals[i];
        path_handle_t source_path_handle = graph->get_path_handle_of_step(interval.first);

        // Skip empty intervals
        if (interval.first == graph->path_end(source_path_handle)) {
            continue;
        }

        // Look up pre-computed source path offsets
        int64_t source_start = step_offset_cache.at(interval.first);
        int64_t source_end = step_offset_cache.at(interval.second);

        // Find reference coordinates using snarl bounds or fallback to BFS
        string base_path_name;
        string ref_path_name = ".";
        int64_t ref_start = 0;
        int64_t ref_end = 0;

        const pair<nid_t, nid_t>& snarl_bounds = this->interval_snarl_bounds[i];

        if (snarl_bounds.first != 0 && snarl_bounds.second != 0 &&
            ref_node_positions.count(snarl_bounds.first) && ref_node_positions.count(snarl_bounds.second)) {
            // Use snarl boundary nodes for reference coordinates
            auto& left_pos = ref_node_positions[snarl_bounds.first];
            auto& right_pos = ref_node_positions[snarl_bounds.second];

            // Both boundary nodes should be on the same reference path
            path_handle_t ref_path_handle = left_pos.first;
            ref_path_name = graph->get_path_name(ref_path_handle);

            int64_t left_node_len = graph->get_length(graph->get_handle(snarl_bounds.first));
            int64_t right_node_len = graph->get_length(graph->get_handle(snarl_bounds.second));

            int64_t left_start = left_pos.second;
            int64_t left_end = left_start + left_node_len;
            int64_t right_start = right_pos.second;
            int64_t right_end = right_start + right_node_len;

            ref_start = min(left_start, right_start);
            ref_end = max(left_end, right_end);

            // Get base path name from the reference path
            base_path_name = ref_path_name;
            subrange_t subrange;
            base_path_name = Paths::strip_subrange(base_path_name, &subrange);
        } else {
            // Fallback to BFS-based get_reference_nodes() for sentinel intervals
            nid_t first_node = graph->get_id(graph->get_handle_of_step(interval.first));
            vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, true);

            if (!ref_nodes.empty()) {
                nid_t ref_node_id = ref_nodes.at(0).second;
                int64_t ref_interval_idx = this->node_to_interval.at(ref_node_id);
                path_handle_t ref_path_handle = graph->get_path_handle_of_step(augref_intervals[ref_interval_idx].first);
                base_path_name = graph->get_path_name(ref_path_handle);
                subrange_t subrange;
                base_path_name = Paths::strip_subrange(base_path_name, &subrange);

                ref_path_name = graph->get_path_name(ref_path_handle);

                // Find the offset of ref_node_id in the reference path
                if (ref_node_positions.count(ref_node_id)) {
                    ref_start = ref_node_positions[ref_node_id].second;
                    ref_end = ref_start + graph->get_length(graph->get_handle(ref_node_id));
                }
            } else {
                // Fallback to source path
                string source_path_name = graph->get_path_name(source_path_handle);
                subrange_t subrange;
                base_path_name = Paths::strip_subrange(source_path_name, &subrange);
            }
        }

        // If augref_sample_name is set, replace sample in base_path_name
        if (!augref_sample_name.empty()) {
            PathSense sense;
            string sample, locus;
            size_t haplotype, phase_block;
            subrange_t subrange;
            PathMetadata::parse_path_name(base_path_name, sense, sample, locus, haplotype, phase_block, subrange);

            if (sample.empty()) {
                base_path_name = augref_sample_name + "#0#" + base_path_name;
            } else {
                base_path_name = PathMetadata::create_path_name(sense, augref_sample_name, locus, haplotype, phase_block, subrange);
            }
        }

        // Get next augref index for this base path
        int64_t augref_index = ++local_augref_counter[base_path_name];
        string augref_name = make_augref_name(base_path_name, augref_index);

        // Resolve source path name to full-path coordinates using cached result.
        // Parses path name once per source path to extract subrange offset and
        // strip the #0 phase block.
        auto cache_it = source_display_cache.find(source_path_handle);
        if (cache_it == source_display_cache.end()) {
            string source_path_name = graph->get_path_name(source_path_handle);
            PathSense sense;
            string sample, locus;
            size_t haplotype, phase_block;
            subrange_t subrange;
            PathMetadata::parse_path_name(source_path_name, sense, sample, locus, haplotype, phase_block, subrange);
            int64_t subrange_offset = (subrange != PathMetadata::NO_SUBRANGE) ? subrange.first : 0;
            if (phase_block == 0) {
                phase_block = PathMetadata::NO_PHASE_BLOCK;
                sense = PathSense::REFERENCE;
            }
            string display_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block, PathMetadata::NO_SUBRANGE);
            cache_it = source_display_cache.emplace(source_path_handle, make_pair(std::move(display_name), subrange_offset)).first;
        }
        const string& display_source_name = cache_it->second.first;
        int64_t display_source_start = source_start + cache_it->second.second;
        int64_t display_source_end = source_end + cache_it->second.second;

        // Output the BED line
        os << display_source_name << "\t"
           << display_source_start << "\t"
           << display_source_end << "\t"
           << augref_name << "\t"
           << ref_path_name << "\t"
           << ref_start << "\t"
           << ref_end << "\n";
    }
}

}
