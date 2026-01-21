#include "altpaths.hpp"
#include <sstream>
#include <algorithm>
#include <queue>
#include <iomanip>
#include <set>

//#define debug

namespace vg {

using namespace std;

const string AltPathsCover::altpath_suffix = "_alt";

string AltPathsCover::make_altpath_name(const string& base_path_name, int64_t alt_index) {
    return base_path_name + altpath_suffix + to_string(alt_index);
}

bool AltPathsCover::is_altpath_name(const string& path_name) {
    // Find the last occurrence of "_alt"
    size_t pos = path_name.rfind(altpath_suffix);
    if (pos == string::npos || pos + altpath_suffix.length() >= path_name.length()) {
        return false;
    }
    // Check that everything after "_alt" is digits
    for (size_t i = pos + altpath_suffix.length(); i < path_name.length(); ++i) {
        if (!isdigit(path_name[i])) {
            return false;
        }
    }
    return true;
}

string AltPathsCover::parse_base_path(const string& altpath_name) {
    if (!is_altpath_name(altpath_name)) {
        return altpath_name;
    }
    size_t pos = altpath_name.rfind(altpath_suffix);
    return altpath_name.substr(0, pos);
}

int64_t AltPathsCover::parse_alt_index(const string& altpath_name) {
    if (!is_altpath_name(altpath_name)) {
        return -1;
    }
    size_t pos = altpath_name.rfind(altpath_suffix);
    return stoll(altpath_name.substr(pos + altpath_suffix.length()));
}

void AltPathsCover::clear(MutablePathMutableHandleGraph* graph) {
    vector<path_handle_t> altpaths_to_remove;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (is_altpath_name(graph->get_path_name(path_handle))) {
            altpaths_to_remove.push_back(path_handle);
        }
    });
    for (path_handle_t path_handle : altpaths_to_remove) {
        graph->destroy_path(path_handle);
    }
}

void AltPathsCover::compute(const PathHandleGraph* graph,
                            SnarlManager* snarl_manager,
                            const unordered_set<path_handle_t>& reference_paths,
                            int64_t minimum_length) {

    // start from scratch
    this->altpath_intervals.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        this->altpath_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                    graph->path_end(ref_path_handle)));
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (node_to_interval.count(node_id)) {
                cerr << "[altpaths error]: node " << node_id << " covered by two reference paths,"
                     << " including " << graph->get_path_name(ref_path_handle) << " and "
                     << graph->get_path_name(graph->get_path_handle_of_step(altpath_intervals.at(node_to_interval.at(node_id)).first))
                     << ". Altpath support currently requires disjoint acyclic reference paths" << endl;
                exit(1);
            }
            node_to_interval[node_id] = altpath_intervals.size() - 1;
        });
    }
    this->num_ref_intervals = this->altpath_intervals.size();

#ifdef debug
#pragma omp critical(cerr)
    cerr << "[altpaths] Selected " << altpath_intervals.size() << " rank=0 reference paths" << endl;
#endif

    // we use the path traversal finder for everything
    PathTraversalFinder path_trav_finder(*graph);

    // we collect the altpath cover in parallel as a list of path fragments
    size_t thread_count = get_thread_count();
    vector<vector<pair<step_handle_t, step_handle_t>>> altpath_intervals_vector(thread_count);
    vector<unordered_map<nid_t, int64_t>> node_to_interval_vector(thread_count);

    // we process top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
        // per-thread output
        vector<pair<step_handle_t, step_handle_t>>& thread_altpath_intervals = altpath_intervals_vector[omp_get_thread_num()];
        unordered_map<nid_t, int64_t>& thread_node_to_interval = node_to_interval_vector[omp_get_thread_num()];

        vector<const Snarl*> queue = {snarl};

        while(!queue.empty()) {
            const Snarl* cur_snarl = queue.back();
            queue.pop_back();

            // get the snarl cover
            compute_snarl(*cur_snarl, path_trav_finder, minimum_length,
                          thread_altpath_intervals,
                          thread_node_to_interval);

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
        cerr << "Adding " << altpath_intervals_vector[t].size() << " altpath intervals from thread " << t << endl;
#endif
        // important to go through function rather than do a raw copy since
        // inter-top-level snarl merging may need to happen
        for (int64_t j = 0; j < altpath_intervals_vector[t].size(); ++j) {
            // the true flag at the end disables the overlap check. since they were computed
            // in separate threads, snarls can overlap by a single node
            const pair<step_handle_t, step_handle_t>& interval = altpath_intervals_vector[t][j];
            if (interval.first != graph->path_end(graph->get_path_handle_of_step(interval.first))) {
                add_interval(this->altpath_intervals, this->node_to_interval, altpath_intervals_vector[t][j], true);
            }
        }
        altpath_intervals_vector[t].clear();
        node_to_interval_vector[t].clear();
    }

    // remove any intervals that were made redundant by add_interval
    defragment_intervals();
}

void AltPathsCover::load(const PathHandleGraph* graph,
                         const unordered_set<path_handle_t>& reference_paths) {
    // start from scratch
    this->altpath_intervals.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (graph->get_is_reverse(graph->get_handle_of_step(step_handle))) {
                cerr << "[altpaths] error: Reversed step " << node_id << " found in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All altpath fragments must be forward-only." << endl;
                exit(1);
            }
            if (node_to_interval.count(node_id)) {
                cerr << "[altpaths] error: Cycle found on node " << node_id << " in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All altpath fragments must be acyclic." << endl;
                exit(1);
            }
            node_to_interval[node_id] = altpath_intervals.size();
        });
        this->altpath_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                    graph->path_end(ref_path_handle)));
    }
    this->num_ref_intervals = this->altpath_intervals.size();

    // load existing altpaths from the graph
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (is_altpath_name(path_name)) {
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                node_to_interval[graph->get_id(graph->get_handle_of_step(step_handle))] = altpath_intervals.size();
            });
            this->altpath_intervals.push_back(make_pair(graph->path_begin(path_handle),
                                                        graph->path_end(path_handle)));
        }
    });
}

void AltPathsCover::apply(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));
#ifdef debug
    cerr << "applying altpath cover with " << this->num_ref_intervals << " ref intervals "
         << " and " << this->altpath_intervals.size() << " total intervals" << endl;
#endif

    // Reset altpath counters for each base path
    base_path_alt_counter.clear();

    // First pass: determine the maximum existing alt index for each base path
    // This ensures we don't overwrite existing altpaths
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_altpath_name(path_name)) {
            string base = parse_base_path(path_name);
            int64_t idx = parse_alt_index(path_name);
            if (base_path_alt_counter.count(base)) {
                base_path_alt_counter[base] = max(base_path_alt_counter[base], idx);
            } else {
                base_path_alt_counter[base] = idx;
            }
        }
    });

    // write the altpaths
    for (int64_t i = this->num_ref_intervals; i < this->altpath_intervals.size(); ++i) {
        // Find the reference path this altpath extends from by tracing back to reference
        nid_t first_node = graph->get_id(graph->get_handle_of_step(altpath_intervals[i].first));
        vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, true);

        // Get the reference path name from the reference node
        string base_path_name;
        if (!ref_nodes.empty()) {
            nid_t ref_node_id = ref_nodes.at(0).second;
            int64_t ref_interval_idx = this->node_to_interval.at(ref_node_id);
            path_handle_t ref_path_handle = graph->get_path_handle_of_step(altpath_intervals[ref_interval_idx].first);
            base_path_name = graph->get_path_name(ref_path_handle);
            // Strip any subrange from the reference path name
            subrange_t subrange;
            base_path_name = Paths::strip_subrange(base_path_name, &subrange);
        } else {
            // Fallback to source path if no reference found (shouldn't happen)
            path_handle_t source_path_handle = mutable_graph->get_path_handle_of_step(altpath_intervals[i].first);
            base_path_name = graph->get_path_name(source_path_handle);
            subrange_t subrange;
            base_path_name = Paths::strip_subrange(base_path_name, &subrange);
        }

        // Get next available alt index for this base path
        int64_t alt_index = ++base_path_alt_counter[base_path_name];

        // Create the altpath name
        string altpath_name = make_altpath_name(base_path_name, alt_index);

        // Create the path as REFERENCE sense
        path_handle_t altpath_handle = mutable_graph->create_path_handle(altpath_name, false);

        for (step_handle_t step_handle = altpath_intervals[i].first; step_handle != altpath_intervals[i].second;
             step_handle = mutable_graph->get_next_step(step_handle)) {
            mutable_graph->append_step(altpath_handle, mutable_graph->get_handle_of_step(step_handle));
        }
    }

    this->forwardize_altpaths(mutable_graph);
}

int64_t AltPathsCover::get_rank(nid_t node_id) const {
    // search back to reference in order to find the rank.
    vector<pair<int64_t, nid_t>> ref_steps = this->get_reference_nodes(node_id, true);
    return ref_steps.at(0).first;
}

step_handle_t AltPathsCover::get_step(nid_t node_id) const {
    if (!node_to_interval.count(node_id)) {
        return step_handle_t();
    }
    const pair<step_handle_t, step_handle_t>& interval = altpath_intervals.at(node_to_interval.at(node_id));
    for (step_handle_t step = interval.first; step != interval.second; step = graph->get_next_step(step)) {
        if (graph->get_id(graph->get_handle_of_step(step)) == node_id) {
            return step;
        }
    }
    return step_handle_t();
}

pair<const pair<step_handle_t, step_handle_t>*,
     const pair<step_handle_t, step_handle_t>*>
AltPathsCover::get_parent_intervals(const pair<step_handle_t, step_handle_t>& interval) const {

    pair<const pair<step_handle_t, step_handle_t>*,
         const pair<step_handle_t, step_handle_t>*> parents = make_pair(nullptr, nullptr);

    // since our decomposition is based on snarl traversals, we know that fragments must
    // overlap their parents on snarl end points (at the very least)
    // therefore we can find parents by scanning along the paths.
    step_handle_t left_parent = graph->get_previous_step(interval.first);
    if (left_parent != graph->path_front_end(graph->get_path_handle_of_step(interval.first))) {
        int64_t interval_idx = this->node_to_interval.at(graph->get_id(graph->get_handle_of_step(left_parent)));
        parents.first = &this->altpath_intervals.at(interval_idx);
    }

    step_handle_t right_parent = graph->get_next_step(interval.second);
    if (right_parent != graph->path_end(graph->get_path_handle_of_step(interval.second))) {
        int64_t interval_idx = node_to_interval.at(graph->get_id(graph->get_handle_of_step(right_parent)));
        parents.second = &this->altpath_intervals.at(interval_idx);
    }
    return parents;
}

const vector<pair<step_handle_t, step_handle_t>>& AltPathsCover::get_intervals() const {
    return this->altpath_intervals;
}

const pair<step_handle_t, step_handle_t>* AltPathsCover::get_interval(nid_t node_id) const {
    if (this->node_to_interval.count(node_id)) {
        return &this->altpath_intervals.at(node_to_interval.at(node_id));
    }
    return nullptr;
}

int64_t AltPathsCover::get_num_ref_intervals() const {
    return this->num_ref_intervals;
}

void AltPathsCover::compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                                  vector<pair<step_handle_t, step_handle_t>>& thread_altpath_intervals,
                                  unordered_map<nid_t, int64_t>& thread_node_to_interval) {

    // start by finding the path traversals through the snarl
    vector<vector<step_handle_t>> travs;
    vector<string> trav_names;
    {
        pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs = path_trav_finder.find_path_traversals(snarl);
        travs.reserve(path_travs.first.size());

        // reduce protobuf usage by going back to vector of steps instead of keeping SnarlTraversals around
        for (int64_t i = 0; i < path_travs.first.size(); ++i) {
            string trav_path_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
            if (is_altpath_name(trav_path_name)) {
                // we ignore existing (off-reference) altpaths
#ifdef debug
                cerr << "Warning : skipping existing altpath traversal " << trav_path_name << endl;
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
        add_interval(thread_altpath_intervals, thread_node_to_interval, new_interval);
    }
}

vector<pair<int64_t, int64_t>> AltPathsCover::get_uncovered_intervals(const vector<step_handle_t>& trav,
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

bool AltPathsCover::add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_altpath_intervals,
                                 unordered_map<nid_t, int64_t>& thread_node_to_interval,
                                 const pair<step_handle_t, step_handle_t>& new_interval,
                                 bool global) {

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
            pair<step_handle_t, step_handle_t>& prev_interval = thread_altpath_intervals[prev_idx];
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
    if (new_interval.second != graph->path_end(graph->get_path_handle_of_step(new_interval.second))) {
        nid_t next_node_id = graph->get_id(graph->get_handle_of_step(new_interval.second));
        if (thread_node_to_interval.count(next_node_id)) {
            int64_t next_idx = thread_node_to_interval[next_node_id];
            pair<step_handle_t, step_handle_t>& next_interval = thread_altpath_intervals[next_idx];
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
                        // decomission next_interval
                        next_interval.first = graph->path_end(next_path);
                        next_interval.second = graph->path_front_end(next_path);
                        deleted_idx = next_idx;
                        // extend the previous interval right
                        thread_altpath_intervals[merged_interval_idx].second = new_interval.second;
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
        merged_interval_idx = thread_altpath_intervals.size();
        thread_altpath_intervals.push_back(new_interval);
    }
    for (step_handle_t step = new_interval.first; step != new_interval.second; step = graph->get_next_step(step)) {
        thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = merged_interval_idx;
    }
    if (deleted_idx >= 0) {
        // move the links to the deleted interval to the new interval
        const pair<step_handle_t, step_handle_t>& deleted_interval = thread_altpath_intervals[deleted_idx];
        for (step_handle_t step = deleted_interval.first; step != deleted_interval.second; step = graph->get_next_step(step)) {
            thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = merged_interval_idx;
        }
    }
    return !merged;
}

void AltPathsCover::defragment_intervals() {
    vector<pair<step_handle_t, step_handle_t>> new_intervals;
    this->node_to_interval.clear();
    for (int64_t i = 0; i < this->altpath_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->altpath_intervals[i];
        path_handle_t path_handle = graph->get_path_handle_of_step(interval.first);
        if (interval.first != graph->path_end(path_handle)) {
            new_intervals.push_back(interval);
        }
    }
    this->altpath_intervals = std::move(new_intervals);
    for (int64_t i = 0; i < this->altpath_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->altpath_intervals[i];
        for (step_handle_t step = interval.first; step != interval.second; step = graph->get_next_step(step)) {
            this->node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = i;
        }
    }
}

int64_t AltPathsCover::get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval) {
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

// Ensure all nodes in altpaths are in forward orientation
void AltPathsCover::forwardize_altpaths(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));

    unordered_map<nid_t, nid_t> id_map;
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_altpath_name(path_name)) {
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
                        cerr << "[altpaths] error: Cycle detected in altpath " << path_name << " at node " << mutable_graph->get_id(handle) << endl;
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
        if (is_altpath_name(path_name)) {
            mutable_graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = mutable_graph->get_handle_of_step(step_handle);
                if (mutable_graph->get_is_reverse(handle)) {
                    cerr << "[altpaths] error: Failed to fowardize node " << mutable_graph->get_id(handle) << " in path " << path_name << endl;
                    exit(1);
                }
            });
        }
    });
}

vector<pair<int64_t, nid_t>> AltPathsCover::get_reference_nodes(nid_t node_id, bool first) const {

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

                const pair<step_handle_t, step_handle_t>& altpath_interval = this->altpath_intervals.at(interval_idx);

                // we've hit the reference, fish out its step and stop searching.
                if (interval_idx < this->num_ref_intervals) {
                    output_reference_nodes.push_back(make_pair(distance, current_id));
                    if (first) {
                        break;
                    }
                    continue;
                }

                // search out of the snarl -- any parent traversals will overlap here
                graph->follow_edges(graph->get_handle_of_step(altpath_interval.first), true, [&](handle_t prev) {
                    queue.push(make_pair(distance + 1, graph->get_id(prev)));
                });
                // hack around gbwtgraph bug (feature?) that does not let you decrement path_end
                path_handle_t path_handle = graph->get_path_handle_of_step(altpath_interval.first);
                step_handle_t last_step;
                if (altpath_interval.second == graph->path_end(path_handle)) {
                    last_step = graph->path_back(path_handle);
                } else {
                    last_step = graph->get_previous_step(altpath_interval.second);
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

    assert(output_reference_nodes.size() > 0 && (!first || (output_reference_nodes.size() == 1)));
    return output_reference_nodes;
}

void AltPathsCover::print_stats(ostream& os) {
    // the header
    os << "#BasePath" << "\t"
       << "AltIndex" << "\t"
       << "Length" << "\t"
       << "NodeStart" << "\t"
       << "NodeEnd" << "\t"
       << "Rank" << "\t"
       << "AvgDepth" << endl;

    for (int64_t i = this->num_ref_intervals; i < this->altpath_intervals.size(); ++i) {
        const pair<step_handle_t, step_handle_t>& interval = this->altpath_intervals[i];
        path_handle_t path_handle = graph->get_path_handle_of_step(interval.first);
        string path_name = graph->get_path_name(path_handle);

        int64_t path_length = 0;
        int64_t tot_depth = 0;
        int64_t tot_steps = 0;
        for (step_handle_t step = interval.first; step != interval.second; step = graph->get_next_step(step)) {
            path_length += graph->get_length(graph->get_handle_of_step(step));
            vector<step_handle_t> steps = graph->steps_of_handle(graph->get_handle_of_step(step));
            tot_depth += steps.size();
            ++tot_steps;
        }

        string base_path = parse_base_path(path_name);
        int64_t alt_index = parse_alt_index(path_name);

        nid_t first_node = graph->get_id(graph->get_handle_of_step(interval.first));
        vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, true);
        int64_t rank = ref_nodes.at(0).first;

        // interval is open ended, so we go back to last node
        step_handle_t last_step;
        if (interval.second == graph->path_end(graph->get_path_handle_of_step(interval.first))) {
            last_step = graph->path_back(graph->get_path_handle_of_step(interval.first));
        } else {
            last_step = graph->get_previous_step(interval.second);
        }
        handle_t last_handle = graph->get_handle_of_step(last_step);

        os << base_path << "\t"
           << alt_index << "\t"
           << path_length << "\t"
           << (graph->get_is_reverse(graph->get_handle_of_step(interval.first)) ? "<" : ">")
           << graph->get_id(graph->get_handle_of_step(interval.first)) << "\t"
           << (graph->get_is_reverse(last_handle) ? "<" : ">")
           << graph->get_id(last_handle) << "\t"
           << rank << "\t"
           << std::fixed << std::setprecision(2) << ((double)tot_depth / tot_steps) << "\n";
    }
}

}
