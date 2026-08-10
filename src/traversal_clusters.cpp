#include "traversal_clusters.hpp"
#include "traversal_finder.hpp"
#include "integrated_snarl_finder.hpp"
#include "snarl_distance_index.hpp"
#include "snarls.hpp"
#include "clip.hpp"
#include "algorithms/dfs.hpp"

//#define debug

namespace vg {

using namespace bdsg;

/// Length-weighted similarity of two snarl alleles.  For two alleles that both carry sequence this
/// is the length-weighted Jaccard coefficient it has always been: |A n B| / |A u B| over multisets
/// of interior handles.
///
/// When one of them is a PURE DELETION -- a traversal straight from snarl start to snarl end, so its
/// interior is empty -- Jaccard has nothing to work with: the intersection is 0 whatever the other
/// allele is, so no threshold could ever merge two deletions.  There, and only there, we measure
/// against the site: sim(deletion, X) = 1 - |X| / max(|site|, |X|), the fraction of the site that
/// BOTH alleles delete.  A 59bp and a 60bp deletion of a 60bp site score 59/60; a deletion and the
/// reference score 0.  Restricting it to that case is the point -- scaling EVERY pair by the site
/// would make two unrelated alleles look similar merely because they sit in a big snarl.
///
/// `site_length` is the total length of the site's reference traversal interior; 0 means "unknown"
/// and falls back to the plain pairwise coefficient for every input.  It is a length rather than the
/// traversal itself because it is invariant across a whole clustering pass, and re-summing it per
/// comparison made every comparison against a pure deletion cost O(site) instead of O(allele).
/// Symmetric; in [0,1] structurally rather than by clamping; exactly 1 iff the two multisets are
/// equal, including two identical pure deletions.
double weighted_traversal_similarity(const PathHandleGraph* graph,
                                     const multiset<handle_t>& target,
                                     const multiset<handle_t>& query,
                                     int64_t site_length) {
    vector<handle_t> intersection_handles;
    std::set_intersection(target.begin(), target.end(),
                          query.begin(), query.end(),
                          std::back_inserter(intersection_handles));
    vector<handle_t> union_handles;
    std::set_union(target.begin(), target.end(),
                   query.begin(), query.end(),
                   std::back_inserter(union_handles));

    auto total_length = [&](const vector<handle_t>& handles) {
        int64_t len = 0;
        for (const handle_t& handle : handles) {
            len += graph->get_length(handle);
        }
        return len;
    };
    int64_t isec_size = total_length(intersection_handles);
    int64_t union_size = total_length(union_handles);

    int64_t denom = union_size;
    if (target.empty() || query.empty()) {
        denom = std::max(site_length, union_size);
    }
    if (denom == 0) {
        // Both alleles and the site are empty, so the two alleles ARE the same allele.  Return 1,
        // not 0/0: a NaN fails both the > and the >= comparisons in cluster_traversals and would
        // make the traversal permanently unclusterable -- the same class of bug as the one this
        // function exists to fix.
        return 1.0;
    }
    return (double)(denom - (union_size - isec_size)) / (double)denom;
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
                                       vector<int>& out_child_snarl_to_trav,
                                       const Traversal* site_ref_trav) {
    
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
        // prune the snarl's boundary handles: every traversal shares them.  A pure deletion (size 2)
        // prunes to the empty set, which is correct -- it carries no sequence.  Keeping its
        // boundaries instead, which is what "trav.size() == 2 ? 0 : 1" used to do, gave it a handle
        // set disjoint from every other allele's BY CONSTRUCTION, so its similarity was structurally
        // 0 and no threshold could merge it with anything.
        multiset<handle_t> sorted_trav;
        for (int64_t j = 1; j + 1 < (int64_t)trav.size(); ++j) {
            sorted_trav.insert(trav[j]);
        }
        sorted_traversals[i] = sorted_trav;
    }
    // The site's reference allele, pruned the same way.  It sets the scale a pure deletion is
    // measured on.  It is not itself clustered, so vg call can supply the VCF reference even though
    // it clusters ALTs only.
    int64_t site_length = 0;
    if (site_ref_trav != nullptr) {
        for (int64_t j = 1; j + 1 < (int64_t)site_ref_trav->size(); ++j) {
            site_length += graph->get_length((*site_ref_trav)[j]);
        }
    }

    for (const int& i : traversal_order) {
        const auto& trav = sorted_traversals[i];
        double max_jaccard = 0;
        int64_t max_cluster_idx = -1;
        for (int64_t j = 0; j < clusters.size(); ++j) {
            const auto& cluster_trav = sorted_traversals[clusters[j][0]];
            double jac = weighted_traversal_similarity(graph, trav, cluster_trav, site_length);
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
    unordered_map<handle_t, vector<int>> handle_to_child;
    for (int64_t i = 0; i < child_snarls.size(); ++i) {
        handle_to_child[child_snarls[i].first].push_back(i);
        handle_to_child[child_snarls[i].second].push_back(i);
    }

    // use the index to find which snarls are fully contained in a given traversal
    // this is a linear scan of the traversal, with above index checked (twice) for each handle
    function<vector<int>(const Traversal&, bool)> get_contained_snarls = [&] (const Traversal& trav, bool fully_contained) {
        map<int, int> fw_count;
        map<int, int> rv_count;
        for (const handle_t& handle : trav) {
            if (handle_to_child.count(handle)) {
                for (int child : handle_to_child[handle]) {                    
                    fw_count[child] += 1;
                }
            }
            handle_t rhandle = graph->flip(handle);
            if (handle_to_child.count(rhandle)) {
                for (int child : handle_to_child[handle]) {
                    rv_count[child] += 1;
                }
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

#ifdef debug
    cerr << "snarl " << graph_interval_to_string(graph, start_handle, end_handle) << endl;
#endif
    
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
        cerr << "trav " << i << ": "
             << "n=" << trav_names[i] << " "
             << "i=" << graph_interval_to_string(graph, graph->get_handle_of_step(path_intervals[i].first),
                                                 graph->get_handle_of_step(path_intervals[i].second))
             << " t=" << traversal_to_string(graph, trav)
             << " a=" << allele <<  " r=" << trav_reversed[i] << endl;
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
            const Traversal& canonical_trav = path_travs[canonical_idx];
            Traversal canonical_trav_flip;
            for (auto i = canonical_trav.rbegin(); i != canonical_trav.rend(); ++i) {
                canonical_trav_flip.push_back(graph->flip(*i));
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
                    const Traversal& replacement_trav = trav_reversed[replace_idx] ? canonical_trav_flip : canonical_trav;
#ifdef debug
                    cerr << "editing interval of size " << path_travs[replace_idx].size()
                         << " from path " << trav_names[replace_idx] << " with canonical interval of size "
                         << replacement_trav.size() << " from path "
                         << trav_names[canonical_idx] << endl
                         << "--interval to replace: "
                         << graph_interval_to_string(graph, graph->get_handle_of_step(interval_to_replace.first),
                                                     graph->get_handle_of_step(interval_to_replace.second)) << endl
                         << "--interval coming in: " << traversal_to_string(graph, replacement_trav) << endl;

#endif
                    assert(graph->get_handle_of_step(interval_to_replace.first) == replacement_trav.front());
                    assert(graph->get_handle_of_step(interval_to_replace.second) == replacement_trav.back());
                    graph->rewrite_segment(interval_to_replace.first, graph->get_next_step(interval_to_replace.second),
                                           replacement_trav);
                }
            }
        }
    }
}

void merge_equivalent_traversals_in_graph(MutablePathHandleGraph* graph, const unordered_set<path_handle_t>& selected_paths,
                                          bool use_snarl_manager) {

    // only consider embedded paths that span snarl
    PathTraversalFinder path_trav_finder(*graph);

    if (use_snarl_manager) {
        // compute the snarls using the old snarl manager
        IntegratedSnarlFinder finder(*graph);
        SnarlManager snarl_manager(std::move(finder.find_snarls_parallel()));

        deque<const Snarl*> queue;
        snarl_manager.for_each_top_level_snarl([&](const Snarl* snarl) {
            queue.push_back(snarl);
        });

        while (!queue.empty()) {
            const Snarl* snarl = queue.front();
            queue.pop_front();
            handle_t start_handle = graph->get_handle(snarl->start().node_id(), snarl->start().backward());
            handle_t end_handle = graph->get_handle(snarl->end().node_id(), snarl->end().backward());
            merge_equivalent_traversals_in_snarl(graph, selected_paths, path_trav_finder, start_handle, end_handle);
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
            for (const Snarl* child : children) {
                queue.push_back(child);
            }
        }
    } else {
        // compute the snarls using the distance index
        // this is what we want to do going forward since it uses the new api, no protobuf etc,
        // but unfortunately it seems way slower on some graphs, hence
        SnarlDistanceIndex distance_index;
        {
            IntegratedSnarlFinder snarl_finder(*graph);
            fill_in_distance_index(&distance_index, graph, &snarl_finder, 0);
        }

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

static bool simplify_snarl_using_traversals(MutablePathMutableHandleGraph* graph, PathTraversalFinder& path_trav_finder,
                                            const handle_t& start_handle, const handle_t& end_handle,
                                            unordered_map<path_handle_t, int64_t>& ref_path_to_rank,
                                            const unordered_map<path_handle_t, int64_t>& selected_ref_paths,
                                            int64_t level,
                                            int64_t max_snarl_length,
                                            double min_jaccard,
                                            int64_t min_fragment_length,
                                            unordered_set<nid_t>& nodes_to_remove,
                                            unordered_set<edge_t>& edges_to_remove) {
    
    // find the path traversals through the snarl
    vector<Traversal> path_travs;
    vector<PathInterval> path_intervals;
    std::tie(path_travs, path_intervals) = path_trav_finder.find_path_traversals(start_handle, end_handle);
    vector<string> trav_names(path_travs.size());
    vector<path_handle_t> trav_paths(path_travs.size());
    vector<int64_t> trav_lengths(path_travs.size(), 0);

#ifdef debug
    cerr << "snarl " << graph_interval_to_string(graph, start_handle, end_handle) << endl;
    cerr << "reference paths";
    for (const auto& rp : ref_paths) {
        cerr << " " << graph->get_path_name(rp);
    }
    cerr << endl;
#endif

    // fill out traversal information (copied from merge_equivalent_traversals_in_snarl() above)
    // and find the reference traversal
    int64_t min_rank = numeric_limits<int64_t>::max();
    int64_t alphabetically_first_ref_trav_idx = -1;
    int64_t max_trav_length = 0;
    for (int64_t i = 0; i < path_travs.size(); ++i) {
        const Traversal& trav = path_travs[i];
        trav_paths[i] = graph->get_path_handle_of_step(path_intervals[i].first);
        trav_names[i] = graph->get_path_name(trav_paths[i]);

#ifdef debug
        cerr << "trav " << i << ": "
             << "n=" << trav_names[i] << " "
             << "i=" << graph_interval_to_string(graph, graph->get_handle_of_step(path_intervals[i].first),
                                                 graph->get_handle_of_step(path_intervals[i].second))
             << " t=" << traversal_to_string(graph, trav, 100) << endl;
#endif

        // note: we are excluding snarl boundaries from length calc here
        for (int64_t j = 1; j < trav.size() - 1; ++j) {
            trav_lengths[i] += graph->get_length(trav[j]);
        }
        max_trav_length = max(max_trav_length, trav_lengths[i]);

        // determine level of path (ie are we seeing it now or was it found at a lower
        // snarl level)
        int64_t path_rank = level;
        if (ref_path_to_rank.count(trav_paths[i])) {
            path_rank = ref_path_to_rank[trav_paths[i]];
        } else {
            ref_path_to_rank[trav_paths[i]] = path_rank;                
        }
        min_rank = min(min_rank, path_rank);
    }

    // choose the reference path
    // the heuristic is used is minimum rank (ie path that spans the top-level snarl)
    // which will be the given reference if available since that comes in with rank 0
    // ties are broken with lexicograph order
    // then choosing the bigger path.
    vector<int64_t> ref_trav_candidates;
    for (int64_t i = 0; i < path_travs.size(); ++i) {
        if (ref_path_to_rank[trav_paths[i]] == min_rank) {
            ref_trav_candidates.push_back(i);
        }
    }
    sort(ref_trav_candidates.begin(), ref_trav_candidates.end(), [&](int64_t i, int64_t j) {
        // A path the user actually named with -P wins outright.  Rank cannot separate them: at a
        // top-level snarl every spanning path is first seen there and so ties with the reference at
        // rank == level == 0, leaving the name comparison below to decide.  Whatever wins becomes
        // ref_trav, which seeds the set of nodes and edges kept, so losing here means the -P
        // reference's own nodes are deleted and its path is left in fragments.
        bool i_selected = selected_ref_paths.count(trav_paths[i]) > 0;
        bool j_selected = selected_ref_paths.count(trav_paths[j]) > 0;
        if (i_selected != j_selected) {
            return i_selected;
        }
        if (trav_names[i] == trav_names[j]) {
            return trav_lengths[i] > trav_lengths[j];
        } else {
            return trav_names[i] < trav_names[j];
        }
    });

    // if there are no reference paths, we bail
    // todo: we could relax this by using the alphabetical path
    if (ref_trav_candidates.empty()) {
#ifdef debug
        cerr << "Level-" << level << " Snarl " << graph_interval_to_string(graph, start_handle, end_handle)
             << " has no reference path in {";
        for (const auto& p : ref_paths) {
            cerr << " " << graph->get_path_name(p);
        }
        cerr << "}" << endl;            
#endif
        return false;
    }

    // find the snarl contents
    // we do this because there's no guarantee in general (ex due to clipping) that every
    // node in the snarl is on a traversal.
    // todo: the old snarl manager had functions for all this (deep_contents). should implement
    // something similar in distance index.
    unordered_set<nid_t> snarl_nodes;
    algorithms::dfs(*graph,
                    [&](const handle_t& h) {
                        snarl_nodes.insert(graph->get_id(h));
                    },
                    [&](const handle_t& h) {},
                    {start_handle},
                    {end_handle, graph->flip(start_handle)});
    unordered_set<edge_t> snarl_edges;
    for (nid_t node_id : snarl_nodes) {
        handle_t handle = graph->get_handle(node_id);
        graph->follow_edges(handle, false, [&](const handle_t& next) {
            if (snarl_nodes.count(graph->get_id(next))) {
                snarl_edges.insert(graph->edge_handle(handle, next));
            }
        });
        graph->follow_edges(handle, true, [&](const handle_t& prev) {
            if (snarl_nodes.count(graph->get_id(prev))) {
                snarl_edges.insert(graph->edge_handle(prev, handle));
            }
        });
    }
                             
    if (path_travs.empty()) {
        cerr << "Snarl " << graph_interval_to_string(graph, start_handle, end_handle) << " contains "
             << snarl_nodes.size() << " nodes, but no path traversals span it and therefore it cannot be simplified" << endl;
        return false;
    }

    // these are the nodes we want to keep
    unordered_set<nid_t> ref_nodes;
    unordered_set<edge_t> ref_edges;
    const Traversal& ref_trav = path_travs[ref_trav_candidates[0]];
    for (int64_t i = 0; i < ref_trav.size(); ++i) {
        ref_nodes.insert(graph->get_id(ref_trav[i]));
        if (i > 0) {
            ref_edges.insert(graph->edge_handle(ref_trav[i-1], ref_trav[i]));
        }
    }

    bool simplify = false;

    // do the length simplification
    if (max_trav_length < max_snarl_length) {
        simplify = true;
    } else if (min_jaccard < 1.0) {
        // do snarl-clustering simplification. note this is only done
        // at the top level.  doing it nested can cause all sorts of weird conflicts.
        set<int64_t> ref_set;
        vector<int> traversal_order;
        for (const int64_t& rc: ref_trav_candidates) {
            traversal_order.push_back((int)rc);
            ref_set.insert(rc);
        }
        for (int i = 0; i < path_travs.size(); ++i) {
            if (!ref_set.count(i)) {
                traversal_order.push_back(i);
            }
        }
        vector<pair<double, int64_t>> info;
        vector<int> child_to_trav;
        // The scale a pure deletion is measured on must be a path the user actually named with -P.
        // ref_trav_candidates[0] is not: at a top-level snarl, level == 0 is the same rank a -P path
        // enters with, so every spanning path ties and the winner is just the alphabetically first
        // name -- which would make the merge depend on how the haplotypes happen to be named.
        // nullptr (no -P path spans this snarl) falls back to pairwise scoring, which merges less.
        const Traversal* site_trav = nullptr;
        for (const int64_t& rc : ref_trav_candidates) {
            if (selected_ref_paths.count(trav_paths[rc])) {
                site_trav = &path_travs[rc];
                break;
            }
        }
        vector<vector<int>> trav_clusters = cluster_traversals(graph, path_travs, traversal_order, {},
                                                               min_jaccard, info, child_to_trav,
                                                               site_trav);
        for (const pair<double, int64_t>& cluster_info : info) {
            if (cluster_info.first < 1.0) {
                simplify = true;
            }
        }

        // some clustering was done, we flag all cluster references to keep, and remove everything else
        assert(trav_clusters[0][0] == ref_trav_candidates[0]); // we've already added ref cluster
        for (int i = 1; i < trav_clusters.size(); ++i) {            
            const Traversal& cluster_ref_trav = path_travs[trav_clusters[i][0]];
            if (simplify) {
                for (int64_t j = 0; j < cluster_ref_trav.size(); ++j) {
                    ref_nodes.insert(graph->get_id(cluster_ref_trav[j]));
                    if (j > 0) {
                        ref_edges.insert(graph->edge_handle(cluster_ref_trav[j-1], cluster_ref_trav[j]));
                    }
                }
            }
        }
    }

    int64_t node_removed_count = 0;
    int64_t edge_removed_count = 0;
    if (simplify) {
        // if simplifyin snarl, flag everything non-reference for removal
        for (const nid_t& node_id : snarl_nodes) {
            if (!ref_nodes.count(node_id)) {
                nodes_to_remove.insert(node_id);
                ++node_removed_count;
            }
        }
        for (const edge_t& edge : snarl_edges) {
            if (!ref_edges.count(edge)) {
                edges_to_remove.insert(edge);
                ++edge_removed_count;
            }
        }
    }
#ifdef debug
    if (node_removed_count + edge_removed_count > 0) {
        cerr << "simplification will remove " << node_removed_count << "  / " << snarl_nodes.size() << " nodes and "
             << edge_removed_count << " edges for snarl "
             << graph_interval_to_string(graph, start_handle, end_handle) << endl;
    }
#endif

    return simplify;
}

void simplify_graph_using_traversals(MutablePathMutableHandleGraph* graph, const string& ref_path_prefix,
                                     int64_t max_snarl_length,
                                     double min_jaccard,
                                     int64_t max_iterations,
                                     int64_t min_fragment_length) {
    // only consider embedded paths that span snarl
    PathTraversalFinder path_trav_finder(*graph);

    // load up the reference paths
    unordered_map<path_handle_t, int64_t> ref_paths;
    std::unordered_map<nid_t, size_t> extra_node_weight;
    constexpr size_t EXTRA_WEIGHT = 10000000000;
    graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC}, [&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (path_name.compare(0, ref_path_prefix.length(), ref_path_prefix) == 0) {
            ref_paths[path_handle] = 0;
            extra_node_weight[graph->get_id(graph->get_handle_of_step(graph->path_begin(path_handle)))] += EXTRA_WEIGHT;
            extra_node_weight[graph->get_id(graph->get_handle_of_step(graph->path_back(path_handle)))] += EXTRA_WEIGHT;
        }
    });

    if (ref_paths.empty()) {
        cerr << "error[vg simplify]: no paths with prefix, " << ref_path_prefix << ", found."
             << " Path simplication presently required at least one reference path to be selected with -P"
             << endl;
        exit(1);
    }

    // both types of normalization selected. we're going to double the iterations
    // and alternate between them
    int64_t input_max_snarl_length = max_snarl_length;
    double input_min_jaccard = min_jaccard;
    int64_t empty_count = 0;
    bool alternate = max_snarl_length > 0 && min_jaccard < 1.0;
    if (alternate) {
        max_iterations *= 2;        
    }
    
    for (int64_t iteration = 0; iteration < max_iterations; ++iteration) {

        // both types of normalization selected. we're going to double the iterations
        // and alternate between them        
        if (alternate) {
            if (iteration % 2 == 0) {
                max_snarl_length = input_max_snarl_length;
                min_jaccard = 1;
            } else {
                max_snarl_length = 0;
                min_jaccard = input_min_jaccard;
            }
        }

        // compute the distance index
        SnarlDistanceIndex distance_index;
        {
            IntegratedSnarlFinder snarl_finder(*graph, extra_node_weight);
            fill_in_distance_index(&distance_index, graph, &snarl_finder, 0);
        }
        
        // do every snarl top-down.
        // to do: there is a potential for overkill - ie a node can get delete by smoothing a top level
        // snarl, then again by a child.  could refactor to be more clever, but I don't yet know
        // if it would save much time in practice. 
        net_handle_t root = distance_index.get_root();
        deque<tuple<net_handle_t, int64_t, unordered_map<path_handle_t, int64_t>>> queue = {make_tuple(root, -1, ref_paths)};
        // we remove all the nodes in one batch to avoid collisions / unececssary path updates
        // todo: does this also need streamlining? also: this setup allows us to work in parallel
        // which could be a possible speedup.    
        unordered_set<nid_t> nodes_to_remove;
        unordered_set<edge_t> edges_to_remove;
    
        while (!queue.empty()) {
            net_handle_t net_handle;
            int64_t level;
            unordered_map<path_handle_t, int64_t> cur_ref_paths;
            std::tie(net_handle, level, cur_ref_paths) = queue.front();
            queue.pop_front();
            bool was_simplified = false;
            if (distance_index.is_snarl(net_handle)) {
                net_handle_t start_bound = distance_index.get_bound(net_handle, false, true);
                net_handle_t end_bound = distance_index.get_bound(net_handle, true, false);
                handle_t start_handle = distance_index.get_handle(start_bound, graph);
                handle_t end_handle = distance_index.get_handle(end_bound, graph);
                was_simplified = simplify_snarl_using_traversals(graph, path_trav_finder, start_handle, end_handle,
                                                                 cur_ref_paths, ref_paths, level,
                                                                 max_snarl_length, min_jaccard,
                                                                 min_fragment_length, nodes_to_remove, edges_to_remove);
            }        
            if (!was_simplified && (
                    net_handle == root || distance_index.is_snarl(net_handle) || distance_index.is_chain(net_handle))) {
                distance_index.for_each_child(net_handle, [&](net_handle_t child_handle) {
                    int64_t next_level = distance_index.is_snarl(child_handle) ? level + 1 : level;
                    queue.push_back(make_tuple(child_handle, next_level, cur_ref_paths));
                });
            }
        }

        if (!nodes_to_remove.empty() || !edges_to_remove.empty()) {
            cerr << "iteration " << iteration;
            if (max_snarl_length > 0) {
                cerr << " (filter snarls < " << max_snarl_length;
            } else {
                cerr << " (merging traversals with similarity > " << min_jaccard;
            }
            cerr << "): deleting " << nodes_to_remove.size() << " nodes and "
                 << edges_to_remove.size() << " edges" << endl;
            // delete the nodes
            delete_nodes_and_chop_paths(graph, nodes_to_remove, edges_to_remove, min_fragment_length);
            empty_count = 0;
        } else {
            ++empty_count;
        }
        if ((alternate && empty_count > 1) || (!alternate && empty_count > 0)) {
            break;
        }
    }
}

unordered_map<string, vector<int>> add_star_traversals(vector<Traversal>& travs,
                                                       vector<string>& names,
                                                       vector<vector<int>>& trav_clusters,
                                                       vector<pair<double, int64_t>>& trav_cluster_info,
                                                       const unordered_map<string, vector<int>>& parent_haplotypes,
                                                       const set<string>& sample_names) {
    // Build a map of what haplotypes are already in the traversals
    unordered_map<string, vector<int>> sample_to_haps;

    assert(names.size() == travs.size());
    for (int64_t i = 0; i < names.size(); ++i) {
        string sample_name = PathMetadata::parse_sample_name(names[i]);
        // for backward compatibility
        if (sample_name.empty()) {
            sample_name = names[i];
        }
        auto phase = PathMetadata::parse_haplotype(names[i]);
        if (!sample_name.empty() && phase == PathMetadata::NO_HAPLOTYPE) {
            // This probably won't fit in an int. Use 0 instead.
            phase = 0;
        }
        sample_to_haps[sample_name].push_back(phase);
    }

    // Find everything that's in parent_haplotypes but not the traversals,
    // and add in dummy star-alleles for them
    for (const auto& parent_sample_haps : parent_haplotypes) {
        string parent_sample_name = PathMetadata::parse_sample_name(parent_sample_haps.first);
        if (parent_sample_name.empty()) {
            parent_sample_name = parent_sample_haps.first;
        }
        if (!sample_names.count(parent_sample_name)) {
            // don't bother for purely reference samples -- we don't need to force an allele for them.
            continue;
        }
        for (int parent_hap : parent_sample_haps.second) {
            bool found = false;
            if (sample_to_haps.count(parent_sample_haps.first)) {
                // note: this is brute-force search, but number of haplotypes usually tiny.
                for (int hap : sample_to_haps[parent_sample_haps.first]) {
                    if (parent_hap == hap) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                travs.push_back(Traversal());
                names.push_back(PathMetadata::create_path_name(PathSense::REFERENCE,
                                                               parent_sample_haps.first,
                                                               "star",
                                                               parent_hap,
                                                               PathMetadata::NO_PHASE_BLOCK,
                                                               PathMetadata::NO_SUBRANGE));
                sample_to_haps[parent_sample_haps.first].push_back(parent_hap);
                trav_clusters.push_back({(int)travs.size() - 1});
                trav_cluster_info.push_back(make_pair(0, 0));
            }
        }
    }

    return sample_to_haps;
}


}
