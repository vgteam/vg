#include "traversal_clusters.hpp"
#include "traversal_finder.hpp"
#include "integrated_snarl_finder.hpp"
#include "snarl_distance_index.hpp"

//#define debug

namespace vg {

using namespace bdsg;

// specialized version of jaccard_coefficient() that weights jaccard by node lenths. 
double weighted_jaccard_coefficient(const PathHandleGraph* graph,
                                    const multiset<handle_t>& target,
                                    const multiset<handle_t>& query) {
    vector<handle_t> intersection_handles;
    std::set_intersection(target.begin(), target.end(),
                          query.begin(), query.end(),
                          std::back_inserter(intersection_handles));
    vector<handle_t> union_handles;
    std::set_union(target.begin(), target.end(),
                   query.begin(), query.end(),
                   std::back_inserter(union_handles));

    int64_t isec_size = 0;
    for (const handle_t& handle : intersection_handles) {
        isec_size += graph->get_length(handle);
    }
    int64_t union_size = 0;
    for (const handle_t& handle : union_handles) {
        union_size += graph->get_length(handle);
    }    
    return (double)isec_size / (double)union_size;
}


vector<int> get_traversal_order(const PathHandleGraph* graph,
                                const vector<Traversal>& traversals,
                                const vector<string>& trav_path_names,
                                const vector<int>& ref_travs,
                                int64_t ref_trav_idx, 
                                const vector<bool>& use_traversal) {
    set<int> ref_set(ref_travs.begin(), ref_travs.end());

    // we want to avoid cyclic references whenever possible, so we downweight in ranking
    vector<int> trav_to_occurrences(traversals.size(), 1);
    unordered_map<string, int> sample_to_occurrence;
    for (int i = 0; i < traversals.size(); ++i) {
        if (use_traversal[i] && !trav_path_names[i].empty()) {
            sample_to_occurrence[PathMetadata::parse_sample_name(trav_path_names[i])] += 1;
        }
    }
    for (int i = 0; i < traversals.size(); ++i) {
        if (use_traversal[i] && !trav_path_names[i].empty()) {
            trav_to_occurrences[i] = sample_to_occurrence[trav_path_names[i]];
        }
    }

    function<bool(int, int)> trav_less = [&](int i, int j) {
        // give reference set preference priority
        bool iref = ref_set.count(i);
        bool jref = ref_set.count(j);
        if (iref != jref) {
            return iref;
        }
        // fall back to occurrences then name comparison 
        if (trav_to_occurrences[i] < trav_to_occurrences[j] ||
            (trav_to_occurrences[i] == trav_to_occurrences[j] && !trav_path_names[i].empty() && (trav_path_names[j].empty() || trav_path_names[i] < trav_path_names[j]))) {
            return true;
        }
        return false;
    };

    vector<int> sorted_travs;
    sorted_travs.reserve(traversals.size());
    for (int64_t i = 0; i < traversals.size(); ++i) {
        if (use_traversal[i] && (i == ref_trav_idx || !ref_set.count(i))) {
            sorted_travs.push_back(i);
        }
    }
#ifdef debug
    cerr << "names:";
    for (auto x : trav_path_names) cerr << " " << x;
    cerr << endl << "ref set:";
    for (auto x : ref_set) cerr << " " << x;
    cerr << endl << "before sort:";
    for (auto x : sorted_travs) cerr << " " << x;
#endif
    std::sort(sorted_travs.begin(), sorted_travs.end(), trav_less);
#ifdef debug
    cerr << endl << "after sort:";
    for (auto x : sorted_travs) cerr << " " << x;
    cerr << endl;
#endif
    assert(ref_travs.empty() || sorted_travs.empty() || ref_set.count(sorted_travs.front()));
    return sorted_travs;
}

vector<vector<int>> cluster_traversals(const PathHandleGraph* graph,
                                       const vector<Traversal>& traversals,
                                       const vector<int>& traversal_order,
                                       const vector<pair<handle_t, handle_t>>& child_snarls,                                       
                                       double min_jaccard,
                                       vector<pair<double, int64_t>>& out_info,
                                       vector<int>& out_child_snarl_to_trav) {
    
    assert(traversal_order.size() <= traversals.size());
    
    // the values are indexes in the input traversals vector
    // the "reference" traversal of each cluster (to which distance is computed)
    // is always its first element
    vector<vector<int>> clusters;

    out_info.resize(traversals.size(), make_pair(-1., 0));
    out_child_snarl_to_trav.resize(child_snarls.size(), -1);

    // keep track of which traversals cover which child snarls
    vector<vector<int>> trav_to_child_snarls = assign_child_snarls_to_traversals(graph, traversals, child_snarls);
    // keep track of which child snarls are covered
    unordered_set<int> uncovered_child_snarls;
    for (const vector<int>& child_snarls : trav_to_child_snarls) {
        for (int i : child_snarls) {
            uncovered_child_snarls.insert(i);
        }
    }
    
    // need the clusters as sorted lists. we'll forget the endpoints while we're at
    // it since they're always shared.  note we work with multisets since we want to
    // count differences between, say, cycle copy numbers. 
    vector<multiset<handle_t>> sorted_traversals(traversals.size());
    for (const int& i : traversal_order) {
        const auto& trav = traversals[i];
        assert(trav.size() >=2);
        // prune snarl nodes as they're always the same
        // todo: does jaccard properly handle empty sets?
        int64_t first = trav.size() == 2 ? 0 : 1;
        int64_t last = trav.size() == 2 ? trav.size() : trav.size() - 1;
        multiset<handle_t> sorted_trav;
        for (int64_t i = first; i < last; ++i) {
            sorted_trav.insert(trav[i]);
        }
        sorted_traversals[i] = sorted_trav;
    } 

    for (const int& i : traversal_order) {
        const auto& trav = sorted_traversals[i];
        double max_jaccard = 0;
        int64_t max_cluster_idx = -1;
        for (int64_t j = 0; j < clusters.size(); ++j) {
            const auto& cluster_trav = sorted_traversals[clusters[j][0]];
            double jac = weighted_jaccard_coefficient(graph, trav, cluster_trav);
            if (jac > max_jaccard) {
                max_jaccard = jac;
                max_cluster_idx = j;
                if (jac == 1) {
                    break;
                }                    
            }        
        }
        if (max_cluster_idx >= 0 && max_jaccard >= min_jaccard) {
            // we've found a suitably similar cluster, add it do that
            clusters[max_cluster_idx].push_back(i);
            out_info[i] = make_pair(max_jaccard, 0);
        } else {
            // there's no cluster close enough, need to start a new one
            clusters.push_back({i});
            out_info[i] = make_pair(1.0, 0);
            // check off all the child snarls it covers
            for (const int& child_snarl_idx : trav_to_child_snarls[i]) {
                if (uncovered_child_snarls.count(child_snarl_idx)) {
                    uncovered_child_snarls.erase(child_snarl_idx);
                    out_child_snarl_to_trav[child_snarl_idx] = i;
                }
            }
        }        
    }

    // break up clusters until all child snarls are covered
    // todo: can/should we factor child coverage into objective of original clustering??
    // todo todo: is it simpler to relax things so that a traversal can be used as a reference
    //            in a nested snarl, even if it doesn't have an exact allele (ie is a cluster reference)
    //            in the parent? -- for now we keep things simple -- all reference alleles are explictly represented in vcf
    vector<vector<int>> new_clusters;
    for (int64_t i = clusters.size() - 1; i >= 0 && !uncovered_child_snarls.empty(); --i) {
        for (int64_t j = clusters[i].size() -1; j > 0 && !uncovered_child_snarls.empty(); --j) {
            const vector<int>& trav_childs = trav_to_child_snarls[clusters[i][j]];
            bool uncovered = false;
            for (int k : trav_childs) {                
                if (uncovered_child_snarls.count(k)) {
                    uncovered = true;
                    uncovered_child_snarls.erase(k);
                    out_child_snarl_to_trav[k] = clusters[i][j];
                }
            }
            if (uncovered) {
                new_clusters.push_back({clusters[i][j]});
                clusters[i].erase(clusters[i].begin() + j);
            }            
        }
    }
    assert(uncovered_child_snarls.empty() || traversal_order.size() < traversals.size());

    // fill in the size deltas
    for (vector<int>& cluster : clusters) {
        // only non-zero for non-empty clusters
        if (cluster.size() > 1) {
            int64_t cluster_ref_length = -1;
            for (int64_t i = 1; i < cluster.size(); ++i) {
                if (out_info[cluster[i]].first < 1) {
                    // get the cluster reference length on-demand
                    if (cluster_ref_length == -1) {
                        cluster_ref_length = 0;
                        for (const handle_t& handle : traversals[cluster[0]]) {
                            cluster_ref_length += graph->get_length(handle);
                        }
                    }
                    // compute the length of the non-ref traversal
                    int64_t length = 0;
                    for (const handle_t& handle : traversals[cluster[i]]) {
                        length += graph->get_length(handle);
                    }
                    // store the delta
                   out_info[cluster[i]].second = length - cluster_ref_length;
                }
                
            }
        }
    }

    return clusters;
}

vector<vector<int>> assign_child_snarls_to_traversals(const PathHandleGraph* graph,
                                                      const vector<Traversal>& traversals,
                                                      const vector<pair<handle_t, handle_t>>& child_snarls) {

    // index the child snarls
    unordered_map<handle_t, int> handle_to_child;
    for (int64_t i = 0; i < child_snarls.size(); ++i) {
        handle_to_child[child_snarls[i].first] = i;
        handle_to_child[child_snarls[i].second] = i;
    }

    // use the index to find which snarls are fully contained in a given traversal
    // this is a linear scan of the traversal, with above index checked (twice) for each handle
    function<vector<int>(const Traversal&, bool)> get_contained_snarls = [&] (const Traversal& trav, bool fully_contained) {
        map<int, int> fw_count;
        map<int, int> rv_count;
        for (const handle_t& handle : trav) {
            if (handle_to_child.count(handle)) {
                fw_count[handle_to_child[handle]] += 1;
            }
            handle_t rhandle = graph->flip(handle);
            if (handle_to_child.count(rhandle)) {
                rv_count[handle_to_child[rhandle]] += 1;
            }
        }
        vector<int> contained_snarls;
        for (const auto& cs_count : fw_count) {
            assert(cs_count.second == 1 || cs_count.second >= 2);
            if (cs_count.second >= 2 || (!fully_contained && cs_count.second == 1)) {
                contained_snarls.push_back(cs_count.first);
            }
        }
        for (const auto& cs_count : rv_count) {
            assert(cs_count.second == 1 || cs_count.second >= 2);
            if (cs_count.second >= 2 || (!fully_contained && cs_count.second == 1)) {            
                contained_snarls.push_back(cs_count.first);
            }
        }        
        return contained_snarls;
    };

    // fill in the output map
    vector<vector<int>> trav_to_child_snarls(traversals.size());

    if (!child_snarls.empty()) {
        for (int64_t i = 0; i < traversals.size(); ++i) {
            trav_to_child_snarls[i] = get_contained_snarls(traversals[i], true);
        }        
    }

    return trav_to_child_snarls;
}

static void merge_equivalent_traversals_in_snarl(MutablePathHandleGraph* graph, const unordered_set<path_handle_t>& selected_paths,
                                                 PathTraversalFinder& path_trav_finder, 
                                                 const handle_t& start_handle, const handle_t& end_handle ) {
    // find the path traversals through the snarl
    vector<Traversal> path_travs;
    vector<PathInterval> path_intervals;
    std::tie(path_travs, path_intervals) = path_trav_finder.find_path_traversals(start_handle, end_handle);
    vector<string> trav_names(path_travs.size());
    vector<path_handle_t> trav_paths(path_travs.size());
    vector<bool> trav_reversed(path_travs.size());

    // organize paths by their alleles strings
    // todo: we could use a fancier index to save memory here (prob not necessary tho)
    unordered_map<string, vector<int64_t>> string_to_travs;
    for (int64_t i = 0; i < path_travs.size(); ++i) {
        const Traversal& trav = path_travs[i];
        string allele;
        for (int64_t j = 1; j < trav.size() - 1; ++j) {
            allele += toUppercase(graph->get_sequence(trav[j]));
        }
        string_to_travs[allele].push_back(i);
        trav_paths[i] = graph->get_path_handle_of_step(path_intervals[i].first);
        trav_names[i] = graph->get_path_name(trav_paths[i]);
        trav_reversed[i] = graph->get_is_reverse(graph->get_handle_of_step(path_intervals[i].first)) !=
            graph->get_is_reverse(start_handle);
#ifdef debug
        cerr << "trav " << i << ": ";
        cerr << "n=" << trav_names[i] << " "
             << "i=" << (graph->get_is_reverse(graph->get_handle_of_step(path_intervals[i].first)) ? "<" : ">")
             << graph->get_id(graph->get_handle_of_step(path_intervals[i].first))
             << (graph->get_is_reverse(graph->get_handle_of_step(path_intervals[i].second)) ? "<" : ">")
             << graph->get_id(graph->get_handle_of_step(path_intervals[i].second)) << " t=";
        for (auto xx : trav) {
            cerr << (graph->get_is_reverse(xx) ? "<" : ">") << graph->get_id(xx);
        }
        cerr << " a=" << allele <<  " r=" << trav_reversed[i] << endl;
#endif
    }        

    // rank the traversals, putting selected ones first otherwise alphabetical
    function<bool(int64_t, int64_t)> trav_idx_less = [&](int64_t i, int64_t j) {
        bool iref = selected_paths.count(trav_paths[i]);
        bool jref = selected_paths.count(trav_paths[j]);
        if (iref != jref) {
            return iref;
        }
        return trav_names[i] < trav_names[j];
    };

    // merge up the paths
    for (const auto& allele_travs : string_to_travs) {
        if (allele_travs.first.length() > 0 && allele_travs.second.size() > 1) {
            const vector<int64_t>& eq_travs = allele_travs.second;
            auto canonical_it = std::min_element(eq_travs.begin(), eq_travs.end(), trav_idx_less);
            assert(canonical_it != eq_travs.end());
            int64_t canonical_idx = *canonical_it;
            Traversal canonical_trav = path_travs[canonical_idx];
            if (trav_reversed[canonical_idx]) {
                Traversal flip_trav;
                for (auto i = canonical_trav.rbegin(); i != canonical_trav.rend(); ++i) {
                    flip_trav.push_back(graph->flip(*i));
                }
                std::swap(canonical_trav, flip_trav);
            }
                
            // edit it into every other path
            //
            // WARNING: in the case of loops, we could be editing the same path more than once
            //          so making the big undocumented assumption that step handles outside edit
            //          are unaffected!!
            for (int64_t i = 0; i < eq_travs.size(); ++i) {
                int64_t replace_idx = eq_travs[i];
                if (replace_idx != canonical_idx && path_travs[replace_idx] != path_travs[canonical_idx]) {
                    PathInterval interval_to_replace = path_intervals[replace_idx];
                    if (trav_reversed[replace_idx]) {
                        std::swap(interval_to_replace.first, interval_to_replace.second);
                    }
#ifdef debug
                    cerr << "editing interval of size " << path_travs[replace_idx].size()
                         << " from path " << trav_names[replace_idx] << " with canonical interval of size "
                         << path_travs[canonical_idx].size() << " from path "
                         << trav_names[canonical_idx] << endl;

                    cerr << "interval to replace first " << graph->get_id(graph->get_handle_of_step(interval_to_replace.first))
                         << " second " << graph->get_id(graph->get_handle_of_step(interval_to_replace.second)) << endl;

#endif                      
                    graph->rewrite_segment(interval_to_replace.first, graph->get_next_step(interval_to_replace.second),
                                           canonical_trav);
                }
            }
        }
    }
}

void merge_equivalent_traversals_in_graph(MutablePathHandleGraph* graph, const unordered_set<path_handle_t>& selected_paths) {  
    // compute the snarls
    SnarlDistanceIndex distance_index;
    {
        IntegratedSnarlFinder snarl_finder(*graph);
        // todo: why can't I pass in 0 below -- I don't want any dinstances!
        fill_in_distance_index(&distance_index, graph, &snarl_finder, 1);
    }

    // only consider embedded paths that span snarl
    PathTraversalFinder path_trav_finder(*graph);

    // do every snarl top-down.  this is because it's possible (tho probably rare) for a child snarl to
    // be redundant after normalizing its parent. don't think the opposite (normalizing child)
    // causes redundant parent.. todo: can we guarantee?!
    net_handle_t root = distance_index.get_root();
    deque<net_handle_t> queue = {root};
    
    while (!queue.empty()) {
        net_handle_t net_handle = queue.front();
        queue.pop_front();
        if (distance_index.is_snarl(net_handle)) {
            net_handle_t start_bound = distance_index.get_bound(net_handle, false, true);
            net_handle_t end_bound = distance_index.get_bound(net_handle, true, false);
            handle_t start_handle = distance_index.get_handle(start_bound, graph);
            handle_t end_handle = distance_index.get_handle(end_bound, graph);
            merge_equivalent_traversals_in_snarl(graph, selected_paths, path_trav_finder, start_handle, end_handle);
        }        
        if (net_handle == root || distance_index.is_snarl(net_handle) || distance_index.is_chain(net_handle)) {
            distance_index.for_each_child(net_handle, [&](net_handle_t child_handle) {
                queue.push_back(child_handle);
            });
        }
    }
}

}
