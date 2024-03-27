#include "rgfa.hpp"
#include <sstream>
#include <algorithm>
#include <queue>

//#define debug

namespace vg {

using namespace std;

const string RGFACover::rgfa_sample_name = "_rGFA_";

string RGFACover::make_rgfa_path_name(const string& path_name, int64_t start, int64_t length,
                                      bool specify_subrange_end) {

    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);
    if (path_subrange == PathMetadata::NO_SUBRANGE) {
        path_subrange = make_pair(0, 0);
    }

    // we move the sample into the locus
    // todo: is there something nicer?
    string rgfa_locus;        
    assert(path_locus != PathMetadata::NO_LOCUS_NAME);    
    if (path_sample != PathMetadata::NO_SAMPLE_NAME) {
        rgfa_locus = path_sample + "::";
    }
    // the contig name will be behind SC...
    rgfa_locus += path_locus;

    // we apply the subrange offset
    path_subrange.first += start;
    path_subrange.second = specify_subrange_end ? path_subrange.first + length : PathMetadata::NO_END_POSITION;


    // and return the final path, with sample/locus/rgfa-rank embedded in locus
    // (as it's a reference path, we alsos strip out the phase block)
    return PathMetadata::create_path_name(PathSense::REFERENCE, RGFACover::rgfa_sample_name,
                                          rgfa_locus, path_haplotype,
                                          PathMetadata::NO_PHASE_BLOCK, path_subrange);

}

bool RGFACover::is_rgfa_path_name(const string& path_name) {
    return PathMetadata::parse_sample_name(path_name) == RGFACover::rgfa_sample_name;
}

pair<string, string> RGFACover::parse_rgfa_locus_name(const string& locus_name) {
    pair<string, string> sample_locus = make_pair(PathMetadata::NO_SAMPLE_NAME, PathMetadata::NO_LOCUS_NAME);

    auto pos = locus_name.find("::");
    if (pos != string::npos) {
        sample_locus.first = locus_name.substr(0, pos);
        sample_locus.second = locus_name.substr(pos + 2); 
    } else {
        sample_locus.second = locus_name;
    }

    return sample_locus;
}    

string RGFACover::revert_rgfa_path_name(const string& rgfa_path_name, bool strip_subrange) {
    if (!is_rgfa_path_name(rgfa_path_name)) {
        return rgfa_path_name;
    }
    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(rgfa_path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);

    std::tie(path_sample, path_locus) = parse_rgfa_locus_name(path_locus);

    if (path_sense == PathSense::REFERENCE && path_sample == PathMetadata::NO_SAMPLE_NAME) {
        path_sense = PathSense::GENERIC;        
    }
    return PathMetadata::create_path_name(path_sense, path_sample,
                                          path_locus, path_haplotype,
                                          path_phase_block, strip_subrange ? PathMetadata::NO_SUBRANGE : path_subrange);
}

void RGFACover::clear(MutablePathMutableHandleGraph* graph) {
    vector<path_handle_t> rgfa_paths;
    graph->for_each_path_of_sample(RGFACover::rgfa_sample_name, [&](path_handle_t path_handle) {
        rgfa_paths.push_back(path_handle);
    });
    for (path_handle_t path_handle : rgfa_paths) {
        graph->destroy_path(path_handle);
    }
}

void RGFACover::compute(const PathHandleGraph* graph,
                        SnarlManager* snarl_manager,
                        const unordered_set<path_handle_t>& reference_paths,
                        int64_t minimum_length) {

    // start from scratch
    this->rgfa_intervals.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        this->rgfa_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                 graph->path_end(ref_path_handle)));
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (node_to_interval.count(node_id)) {
                cerr << "[rgfa error]: node " << node_id << " covered by two reference paths,"
                     << " including " << graph->get_path_name(ref_path_handle)
                     << ". rGFA support current requires disjoint acyclic reference paths" << endl;
                exit(1);
            }
            node_to_interval[node_id] = rgfa_intervals.size() - 1;
        });
    }
    this->num_ref_intervals = this->rgfa_intervals.size();

#ifdef debug
    cerr << "[rgfa] Selected " << rgfa_intervals.size() << " rank=0 reference paths" << endl;
#endif
    
    // we use the path traversal finder for everything
    // (even gbwt haplotypes, as we're using the path handle interface)
    PathTraversalFinder path_trav_finder(*graph, *snarl_manager);

    // we collect the rgfa cover in parallel as a list of path fragments
    size_t thread_count = get_thread_count();
    vector<vector<pair<step_handle_t, step_handle_t>>> rgfa_intervals_vector(thread_count);    
    vector<unordered_map<nid_t, int64_t>> node_to_interval_vector(thread_count);
    
    // we process top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
        // per-thread output
        vector<pair<step_handle_t, step_handle_t>>& thread_rgfa_intervals = rgfa_intervals_vector[omp_get_thread_num()];    
        unordered_map<nid_t, int64_t>& thread_node_to_interval = node_to_interval_vector[omp_get_thread_num()];

        vector<const Snarl*> queue = {snarl};

        cerr << "top level snarl " << pb2json(*snarl) << endl;

        while(!queue.empty()) {
            const Snarl* cur_snarl = queue.back();
            queue.pop_back();

            // get the snarl cover
            compute_snarl(*cur_snarl, path_trav_finder, minimum_length,
                          thread_rgfa_intervals,
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
        int64_t offset = this->rgfa_intervals.size();

        // todo: figure out one-shot stl invocation to move
        for (int64_t j = 0; j < rgfa_intervals_vector[t].size(); ++j) {
            this->rgfa_intervals.push_back(rgfa_intervals_vector[t][j]);
        }
#ifdef debug
        cerr << "Adding " << rgfa_intervals_vector[t].size() << " rgfa intervals from thread " << t << endl;
#endif
        rgfa_intervals_vector[t].clear();

        for (const auto& node_interval : node_to_interval_vector[t]) {
            this->node_to_interval[node_interval.first] = node_interval.second + offset;
        }
        node_to_interval_vector[t].clear();
    }

    // todo: we could optionally go through uncovered nodes and try to greedily search
    // for rgfa intervals at this point, since it seems for some graphs there are
    // regions that don't get found via the traversal finder

}

void RGFACover::load(const PathHandleGraph* graph,
                     const unordered_set<path_handle_t>& reference_paths,
                     bool use_original_paths) {
    // start from scratch
    this->rgfa_intervals.clear();
    this->node_to_interval.clear();
    this->graph = graph;

    // start with the reference paths
    for (const path_handle_t& ref_path_handle : reference_paths) {
        graph->for_each_step_in_path(ref_path_handle, [&](step_handle_t step_handle) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(step_handle));
            if (graph->get_is_reverse(graph->get_handle_of_step(step_handle))) {
                cerr << "[rgfa] error: Reversed step " << node_id << " found in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All rGFA path fragments must be forward-only." << endl;
                exit(1);
            }            
            if (node_to_interval.count(node_id)) {
                cerr << "[rgfa] error: Cycle found on node " << node_id << "in rank-0 reference "
                     << graph->get_path_name(ref_path_handle) << ". All rGFA path fragments must be acyclic." << endl;
                exit(1);

            }
            node_to_interval[node_id] = rgfa_intervals.size();
        });
        this->rgfa_intervals.push_back(make_pair(graph->path_begin(ref_path_handle),
                                                 graph->path_end(ref_path_handle)));
    }
    this->num_ref_intervals = this->rgfa_intervals.size();

    if (!use_original_paths) {
        // just suck up the rgfa fragments right into the data structure
        graph->for_each_path_of_sample(RGFACover::rgfa_sample_name, [&](path_handle_t path_handle) {
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                node_to_interval[graph->get_id(graph->get_handle_of_step(step_handle))] = rgfa_intervals.size();
            });
            this->rgfa_intervals.push_back(make_pair(graph->path_begin(path_handle),
                                                     graph->path_end(path_handle)));
        });
    } else {
        // then the rgfa cover paths
        // since we want to keep our structures in  terms of original paths, we have to map back
        // to them here (if we don't have original paths, then we can't find the overlaps and
        // therefore nesting relationships between them.

        // we start by making a little index in two scans, as I'm worried about quadratic path scan below otherwise
        // (this does not make guarantees about degenerate fragmentation related cases tho)
        unordered_map<pair<string, string>, vector<path_handle_t>> sample_locus_to_paths;    
        graph->for_each_path_of_sample(RGFACover::rgfa_sample_name, [&](path_handle_t path_handle) {
            string locus_name = graph->get_locus_name(path_handle);
            sample_locus_to_paths[parse_rgfa_locus_name(locus_name)] = {};
        });
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            pair<string, string> sample_locus = make_pair(graph->get_sample_name(path_handle), graph->get_locus_name(path_handle));
            if (sample_locus_to_paths.count(sample_locus)) {
                sample_locus_to_paths[sample_locus].push_back(path_handle);
            }
        });

        // next, we scan each rgfa path fragment, and use the index to semi-quickly find its source path interval
        // todo: An inconsistency between cover paths and source paths is a possibility if someone messed up their graph
        // so should probably have a better error message than the asserts below (ie if exact interval match not found)
        graph->for_each_path_of_sample(RGFACover::rgfa_sample_name, [&](path_handle_t path_handle) {
            // pase the rgfa locus to get the original sample and locus
            pair<string, string> source_sample_locus = parse_rgfa_locus_name(graph->get_locus_name(path_handle));
            // find the sample in our index
            const vector<path_handle_t>& source_paths = sample_locus_to_paths.at(source_sample_locus);
            // find the containing path
            subrange_t rgfa_subrange = graph->get_subrange(path_handle);
            assert(rgfa_subrange != PathMetadata::NO_SUBRANGE);
            int64_t rgfa_haplotype = graph->get_haplotype(path_handle);
            // allow match between 0 and NO_HAPLOTYPE
            if (rgfa_haplotype == PathMetadata::NO_HAPLOTYPE) {
                rgfa_haplotype = 0;
            }
            const path_handle_t* source_path = nullptr;
            subrange_t source_subrange;
            for (const path_handle_t& source_path_candidate : source_paths) {
                int64_t source_haplotype = graph->get_haplotype(source_path_candidate);
                if (source_haplotype == PathMetadata::NO_HAPLOTYPE) {
                    source_haplotype = 0;
                }
                if (source_haplotype == rgfa_haplotype) {
                    source_subrange = graph->get_subrange(source_path_candidate);
                    if (source_subrange == PathMetadata::NO_SUBRANGE) {
                        source_subrange.first = 0;
                    }
                    if (source_subrange == PathMetadata::NO_SUBRANGE || source_subrange.second == PathMetadata::NO_END_POSITION) {
                        source_subrange.second = 0;
                        graph->for_each_step_in_path(source_path_candidate, [&](step_handle_t step) {
                            source_subrange.second += graph->get_length(graph->get_handle_of_step(step));
                        });
                    }
                    if (rgfa_subrange.first >= source_subrange.first && rgfa_subrange.second <= source_subrange.second) {
                        source_path = &source_path_candidate;
                        break;
                    }
                }
            }
            assert(source_path != nullptr);
            // now find the exact interval in the containing path and update our data structure
            bool found_start = false;
            step_handle_t source_start;
            int64_t cur_offset = 0;
            graph->for_each_step_in_path(*source_path, [&](step_handle_t cur_step) {
                if (cur_offset + source_subrange.first == rgfa_subrange.first) {
                    source_start = cur_step;
                    found_start = true;
                } else {
                    cur_offset += graph->get_length(graph->get_handle_of_step(cur_step));
                }
                return !found_start;
            });
            assert(found_start);
            assert(graph->get_id(graph->get_handle_of_step(source_start)) ==
                   graph->get_id(graph->get_handle_of_step(graph->path_begin(path_handle))));
            assert(graph->get_is_reverse(graph->get_handle_of_step(source_start)) ==
                   graph->get_is_reverse(graph->get_handle_of_step(graph->path_begin(path_handle))));

            bool found_end = false;
            step_handle_t source_end;
            int64_t num_steps = 0;
            for (step_handle_t cur_step = source_start;; cur_step = graph->get_next_step(cur_step)) {
                ++num_steps;
                cur_offset += graph->get_length(graph->get_handle_of_step(cur_step));
            
                if (cur_offset + source_subrange.first == rgfa_subrange.second) {
                    found_end = true;
                    source_end = cur_step;
                    break;
                }
                if (!graph->has_next_step(cur_step)) {
                    break;
                }
            }
            assert(found_end);
            assert(num_steps == graph->get_step_count(path_handle));

            // we can finally add our interval
            source_end = graph->get_next_step(source_end);
            this->rgfa_intervals.push_back(make_pair(source_start, source_end));
            for (step_handle_t cur_step = source_start; cur_step != source_end; cur_step = graph->get_next_step(cur_step)) {
                int64_t cur_id = graph->get_id(graph->get_handle_of_step(cur_step));
                assert(!node_to_interval.count(cur_id));
                node_to_interval[cur_id] = rgfa_intervals.size() - 1;
            }
        });
    }
}

void RGFACover::apply(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));

    // index paths isued in rgfa cover
    unordered_map<step_handle_t, int64_t> step_to_offset;
    unordered_set<path_handle_t> source_path_set;
    for (int64_t i = this->num_ref_intervals; i < this->rgfa_intervals.size(); ++i) {
        path_handle_t source_path_handle = mutable_graph->get_path_handle_of_step(rgfa_intervals[i].first);
        step_to_offset[rgfa_intervals[i].first] = -1;
        source_path_set.insert(source_path_handle);
    }
    vector<path_handle_t> source_paths(source_path_set.begin(), source_path_set.end()); // for old-style omp interface
#pragma omp parallel for
    for (int64_t i = 0; i < source_paths.size(); ++i) {
        int64_t offset = 0;
        graph->for_each_step_in_path(source_paths[i], [&](step_handle_t step_handle) {
            if (step_to_offset.count(step_handle)) {
                step_to_offset[step_handle] = offset;
            }
            offset += graph->get_length(graph->get_handle_of_step(step_handle));
        });
    }
    
    // compute the offsets in parallel
    vector<int64_t> rgfa_offsets(this->rgfa_intervals.size());
    vector<int64_t> rgfa_lengths(this->rgfa_intervals.size());

#pragma omp parallel for
    for (int64_t i = this->num_ref_intervals; i < this->rgfa_intervals.size(); ++i) {
        rgfa_offsets[i] = step_to_offset.at(rgfa_intervals[i].first);
        rgfa_lengths[i] = 0;
        for (step_handle_t step_handle = rgfa_intervals[i].first; step_handle != rgfa_intervals[i].second;
             step_handle = mutable_graph->get_next_step(step_handle)) {
            rgfa_lengths[i] += graph->get_length(graph->get_handle_of_step(step_handle));
        }        
    }
    step_to_offset.clear();

    // write the rgfa paths
    for (int64_t i = this->num_ref_intervals; i < this->rgfa_intervals.size(); ++i) {
        path_handle_t source_path_handle = mutable_graph->get_path_handle_of_step(rgfa_intervals[i].first);
        string source_path_name = graph->get_path_name(source_path_handle);
        string rgfa_path_name = make_rgfa_path_name(source_path_name, rgfa_offsets[i], rgfa_lengths[i]);
        path_handle_t rgfa_path_handle = mutable_graph->create_path_handle(rgfa_path_name);
        for (step_handle_t step_handle = rgfa_intervals[i].first; step_handle != rgfa_intervals[i].second;
             step_handle = mutable_graph->get_next_step(step_handle)) {
            mutable_graph->append_step(rgfa_path_handle, mutable_graph->get_handle_of_step(step_handle));
        }
    }

    this->forwardize_rgfa_paths(mutable_graph);
}

int64_t RGFACover::get_rank(nid_t node_id) const {

    // search back to reference in order to find the rank.
    vector<pair<int64_t, nid_t>> ref_steps = this->get_reference_nodes(node_id, true);
    return ref_steps.at(0).first;
}
    
pair<const pair<step_handle_t, step_handle_t>*,
     const pair<step_handle_t, step_handle_t>*>
RGFACover::get_parent_intervals(const pair<step_handle_t, step_handle_t>& interval) const {

    pair<const pair<step_handle_t, step_handle_t>*,
         const pair<step_handle_t, step_handle_t>*> parents = make_pair(nullptr, nullptr);
    
    // since our decomposition is baseds on snarl tranversals, we know that fragments must
    // overlap their parents on snarl end points (at the very least)
    // therefore we can find parents by scanning along the rgfa paths.    
    step_handle_t left_parent = graph->get_previous_step(interval.first);
    if (left_parent != graph->path_front_end(graph->get_path_handle_of_step(interval.first))) {
        int64_t interval_idx = this->node_to_interval.at(graph->get_id(graph->get_handle_of_step(left_parent)));
        parents.first = &this->rgfa_intervals.at(interval_idx);
    }

    step_handle_t right_parent = graph->get_next_step(interval.second);
    if (right_parent != graph->path_end(graph->get_path_handle_of_step(interval.second))) {
        int64_t interval_idx = node_to_interval.at(graph->get_id(graph->get_handle_of_step(right_parent)));
        parents.second = &this->rgfa_intervals.at(interval_idx);
    }
    return parents;
}

const vector<pair<step_handle_t, step_handle_t>>& RGFACover::get_intervals() const {
    return this->rgfa_intervals;
}

const pair<step_handle_t, step_handle_t>* RGFACover::get_interval(nid_t node_id) const {
    if (this->node_to_interval.count(node_id)) {
        return &this->rgfa_intervals.at(node_to_interval.at(node_id));
    }
    return nullptr;
}


void RGFACover::compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                              vector<pair<step_handle_t, step_handle_t>>& thread_rgfa_intervals,
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
            if (is_rgfa_path_name(trav_path_name)) {
                // we ignore existing (off-reference) rGFA paths
                // todo: shoulgd there be better error handling?                
                cerr << "Warning : skipping existing rgfa traversal " << trav_path_name << endl;
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
        cerr << "best trav: ";
        for (auto xx : trav) cerr << " " << graph->get_id(graph->get_handle_of_step(xx));
        cerr << endl << "uncovered interval [" << uncovered_interval.first << "," << uncovered_interval.second << "]" <<endl;
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
        cerr << "adding interval with length " << interval_length << endl;
#endif
        add_interval(thread_rgfa_intervals, thread_node_to_interval, new_interval);
    }
}

vector<pair<int64_t, int64_t>> RGFACover::get_uncovered_intervals(const vector<step_handle_t>& trav,
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

bool RGFACover::add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_rgfa_intervals,
                             unordered_map<nid_t, int64_t>& thread_node_to_interval,
                             const pair<step_handle_t, step_handle_t>& new_interval) {

    bool merged = false;
    path_handle_t path_handle = graph->get_path_handle_of_step(new_interval.first);

    // check the before-first step.  if it's in an interval then it must be immediately
    // preceeding so we merge the new interval to the end of the found interval
    step_handle_t before_first_step = graph->get_previous_step(new_interval.first);
    if (before_first_step != graph->path_front_end(graph->get_path_handle_of_step(before_first_step))) {
        nid_t prev_node_id = graph->get_id(graph->get_handle_of_step(before_first_step));
        if (thread_node_to_interval.count(prev_node_id)) {
            pair<step_handle_t, step_handle_t>& prev_interval = thread_rgfa_intervals[thread_node_to_interval[prev_node_id]];
            if (graph->get_path_handle_of_step(prev_interval.first) == path_handle) {
                assert(prev_interval.second == new_interval.first);
                prev_interval.second = new_interval.second;
                merged = true;
            }
        }
    }
    
    // check the end step.  if it's in an interval then it must be immediately
    // following  we merge the new interval to the front of the found interval
    if (new_interval.second != graph->path_end(graph->get_path_handle_of_step(new_interval.second))) {
        nid_t next_node_id = graph->get_id(graph->get_handle_of_step(new_interval.second));
        if (thread_node_to_interval.count(next_node_id)) {
            pair<step_handle_t, step_handle_t>& next_interval = thread_rgfa_intervals[thread_node_to_interval[next_node_id]];
            path_handle_t next_path = graph->get_path_handle_of_step(next_interval.first);
            if (graph->get_path_handle_of_step(next_interval.first) == path_handle) {
                assert(next_interval.first == new_interval.second);
                next_interval.first = new_interval.first;
                merged = true;
            }
        }
    }

    // add the interval to the local (thread safe) structures
    if (!merged) {
        for (step_handle_t step = new_interval.first; step != new_interval.second; step = graph->get_next_step(step)) {
            thread_node_to_interval[graph->get_id(graph->get_handle_of_step(step))] = thread_rgfa_intervals.size();
        }
        thread_rgfa_intervals.push_back(new_interval);
    }
    
    return !merged;
}

int64_t RGFACover::get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval) {
    path_handle_t path_handle = graph->get_path_handle_of_step(trav.front());
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



// copied pretty much verbatem from
// https://github.com/ComparativeGenomicsToolkit/hal2vg/blob/v1.1.2/clip-vg.cpp#L809-L880
void RGFACover::forwardize_rgfa_paths(MutablePathMutableHandleGraph* mutable_graph) {
    assert(this->graph == static_cast<PathHandleGraph*>(mutable_graph));
    
    unordered_map<nid_t, nid_t> id_map;
    mutable_graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = mutable_graph->get_path_name(path_handle);
        if (is_rgfa_path_name(path_name)) {
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
                        cerr << "[rGFA] error: Cycle detected in rGFA path " << path_name << " at node " << mutable_graph->get_id(handle) << endl;
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
        if (is_rgfa_path_name(path_name)) {
            mutable_graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = mutable_graph->get_handle_of_step(step_handle);
                if (mutable_graph->get_is_reverse(handle)) {
                    cerr << "[rGFA] error: Failed to fowardize node " << mutable_graph->get_id(handle) << " in path " << path_name << endl;
                    exit(1);
                }
            });
        }
    });
}

vector<pair<int64_t, nid_t>> RGFACover::get_reference_nodes(nid_t node_id, bool first) const {

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

                const pair<step_handle_t, step_handle_t>& rgfa_interval = this->rgfa_intervals.at(interval_idx);

                // we've hit the reference, fish out its step and stop searching.
                if (interval_idx < this->num_ref_intervals) {
                    output_reference_nodes.push_back(make_pair(distance, current_id));
                    if (first) {
                        break;
                    }
                    continue;
                }

                // search out of the snarl -- any parent traversals will overlap here                
                graph->follow_edges(graph->get_handle_of_step(rgfa_interval.first), true, [&](handle_t prev) {
                    queue.push(make_pair(distance + 1, graph->get_id(prev)));
                });
                // hack around gbwtgraph bug (feature?) that does not let you decrement path_end
                path_handle_t path_handle = graph->get_path_handle_of_step(rgfa_interval.first);
                step_handle_t last_step;
                if (rgfa_interval.second == graph->path_end(path_handle)) {
                    last_step = graph->path_back(path_handle);
                } else {
                    last_step = graph->get_previous_step(rgfa_interval.second);
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

void RGFACover::annotate_vcf(vcflib::VariantCallFile& vcf, ostream& os) {
    vcflib::Variant var(vcf);

    // want a reference position lookup
    unordered_map<nid_t, int64_t> node_to_ref_pos;
    for (int64_t i = 0; i < this->num_ref_intervals; ++i) {
        const pair<step_handle_t, step_handle_t>& ref_interval = this->rgfa_intervals.at(i);
        // assumption: ref intervals span entire path
        assert(graph->path_begin(graph->get_path_handle_of_step(ref_interval.first)) == ref_interval.first);
        int64_t pos = 0;
        for (step_handle_t step_handle = ref_interval.first; step_handle != ref_interval.second;
             step_handle = graph->get_next_step(step_handle)) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            nid_t node_id = graph->get_id(handle);
            assert(!graph->get_is_reverse(handle));            
            assert(!node_to_ref_pos.count(node_id));
            node_to_ref_pos[node_id] = pos;
            pos += graph->get_length(handle);
        }
    }
    
    // also want an rGFA path lookup
    unordered_map<string, path_handle_t> name_to_rgfa_path;
    graph->for_each_path_of_sample(RGFACover::rgfa_sample_name, [&](path_handle_t path_handle) {
        name_to_rgfa_path[RGFACover::revert_rgfa_path_name(graph->get_path_name(path_handle))] = path_handle;
    });

    // remove header lines we're going to add
    vcf.removeInfoHeaderLine("R_CHROM");
    vcf.removeInfoHeaderLine("R_START");
    vcf.removeInfoHeaderLine("R_END");
    vcf.removeInfoHeaderLine("RANK");

    // and add them back, so as not to duplicate them if they are already there
    vcf.addHeaderLine("##INFO=<ID=R_CHROM,Number=1,Type=String,Description=\"Reference Chromsome Name\">");
    vcf.addHeaderLine("##INFO=<ID=R_START,Number=1,Type=String,Description=\"Chromsome Start (0-based inclusive)\">");
    vcf.addHeaderLine("##INFO=<ID=R_END,Number=1,Type=String,Description=\"Reference End (0-based exclusive)\">");
    vcf.addHeaderLine("##INFO=<ID=RANK,Number=1,Type=String,Description=\"rGFA Rank\">");

    os << vcf.header << endl;

    string prev_sequence_name;
    string prev_r_chrom;
    string prev_r_start;
    string prev_r_end;
    string prev_rank;
    while (vcf.getNextVariant(var)) {
        if (name_to_rgfa_path.count(var.sequenceName)) {
            string r_chrom;
            string r_start;
            string r_end;
            string rank;
            if (var.sequenceName == prev_sequence_name) {
                // just use the previous values, which will be the same
                r_chrom = prev_r_chrom;
                r_start = prev_r_start;
                r_end = prev_r_end;
                rank = prev_rank;
            } else {
                // compute from the cover
                path_handle_t path_handle = name_to_rgfa_path.at(var.sequenceName);
                nid_t first_node = graph->get_id(graph->get_handle_of_step(graph->path_begin(path_handle)));
                vector<pair<int64_t, nid_t>> ref_nodes = this->get_reference_nodes(first_node, false);

                int64_t min_ref_pos = numeric_limits<int64_t>::max();
                int64_t max_ref_pos = -1;
                int64_t min_rank = numeric_limits<int64_t>::max();
                for (const pair<int64_t, nid_t>& rank_node : ref_nodes) {
                    step_handle_t ref_step = this->rgfa_intervals.at(node_to_interval.at(rank_node.second)).first;
                    path_handle_t ref_path = graph->get_path_handle_of_step(ref_step);
                    string name = graph->get_path_name(ref_path);
                    // we assume one reference contig (which is built into the whole structure)
                    assert(r_chrom.empty() || r_chrom == name);
                    r_chrom = name;
                    int64_t ref_pos = node_to_ref_pos.at(rank_node.second);
                    // assume snarl is forward on both reference nodes
                    // todo: this won't be exact for some inversion cases, I don't think --
                    //       need to test these and either add check / or move to oriented search
                    min_ref_pos = min(min_ref_pos, ref_pos + (int64_t)graph->get_length(graph->get_handle(rank_node.second)));
                    max_ref_pos = max(max_ref_pos, ref_pos);
                    min_rank = min(min_rank, rank_node.first);
                }

                r_start = std::to_string(min_ref_pos);
                r_end = std::to_string(max_ref_pos);
                rank = std::to_string(min_rank);
            }

            if (!var.info.count("R_CHROM")) {
                var.format.push_back("R_CHROM");
            }
            var.info["R_CHROM"] = {r_chrom};
            if (!var.info.count("R_START")) {
                var.format.push_back("R_START");
            }
            var.info["R_START"] = {r_start};
            if (!var.info.count("R_END")) {
                var.format.push_back("R_END");
            }
            var.info["R_END"] = {r_end};
            if (!var.info.count("RANK")) {
                var.format.push_back("RANK");
            }
            var.info["RANK"] = {rank};
            
            prev_sequence_name = var.sequenceName;
            prev_r_chrom = r_chrom;
            prev_r_start = r_start;
            prev_r_end = r_end;
            prev_rank = rank;
            
        } else {
            //cerr << "rGFA [warning]: VCF reference " << var.sequenceName << " not found in graph" << endl;
        }
        os << var << endl;        
    }
}

}

