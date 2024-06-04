#include "traversal_clusters.hpp"

//#define debug

namespace vg {


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
    
    function<bool(int, int)> trav_less = [&](int i, int j) {
        // give reference set preference priority
        bool iref = ref_set.count(i);
        bool jref = ref_set.count(j);
        if (iref != jref) {
            return iref;
        }
        // fall back to name comparison 
        if (!trav_path_names[i].empty() && (trav_path_names[j].empty() || trav_path_names[i] < trav_path_names[j])) {
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



   


}
