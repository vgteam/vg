#include "alignment_path_offsets.hpp"

//#define debug_mpaln_offsets

namespace vg {
namespace algorithms {

unordered_map<path_handle_t, vector<pair<size_t, bool> > >
alignment_path_offsets(const PathPositionHandleGraph& graph,
                       const Alignment& aln,
                       bool just_min,
                       bool nearby,
                       int64_t search_limit,
                       const std::function<bool(const path_handle_t&)>* path_filter) {
    if (nearby && search_limit == 0) {
        // Fill in the search limit
        search_limit = aln.sequence().size();
    }
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > offsets;
    if (graph.get_path_count() == 0) return offsets;
    for (auto& mapping : aln.path().mapping()) {
        // How many bases does this Mapping cover over?
        size_t mapping_width = mapping_from_length(mapping);
        if (mapping_width == 0 && !nearby) {
            // Just skip over this mapping; it touches no bases.
            continue;
        }
        // We may have to consider both the starts and ends of mappings
        vector<bool> end = {false};
        if (just_min && !nearby) {
            // We want the min actually touched position along each path. It
            // could come from the Mapping start or the Mapping end.
            end.push_back(true);
        }
        // Find the position of this end of this mapping
        pos_t mapping_pos = make_pos_t(mapping.position());
        // Find the positions for this end of this Mapping
        auto pos_offs = algorithms::nearest_offsets_in_paths(&graph, mapping_pos, nearby ? search_limit : -1, path_filter);
        for (auto look_at_end : end) {
            // For the start and the end of the Mapping, as needed
            for (auto& p : pos_offs) {
                // For each path, splice the list of path positions for this Mapping
                // onto the end of the list of positions we found in that path
                auto& v = offsets[p.first];
                for (pair<size_t, bool>& y : p.second) {
                    v.emplace_back(y.second ? y.first - mapping_width : y.first,
                                   y.second);
                }
            }
        }
    }
    if (!nearby && offsets.empty() && search_limit != -1) {
        // find the nearest if we couldn't find any before but we could do a search
        return alignment_path_offsets(graph, aln, just_min, true, search_limit, path_filter);
    }
    if (just_min) {
        // We need the minimum position for each path
        for (auto& p : offsets) {
            auto& v = p.second;
            auto m = *min_element(v.begin(), v.end(),
                                  [](const pair<size_t, bool>& a,
                                     const pair<size_t, bool>& b)
                                  { return a.first < b.first; });
            v.clear();
            v.push_back(m);
        }
    }
    return offsets;
}

unordered_map<path_handle_t, vector<pair<size_t, bool> > >
multipath_alignment_path_offsets(const PathPositionHandleGraph& graph,
                                 const multipath_alignment_t& mp_aln,
                                 const std::function<bool(const path_handle_t&)>* path_filter) {
    
    using path_positions_t = unordered_map<path_handle_t, vector<pair<size_t, bool>>>;
    
    // collect the search results for each mapping on each subpath
    vector<vector<path_positions_t>> search_results(mp_aln.subpath_size());
    for (size_t i = 0; i < mp_aln.subpath_size(); ++i) {
        const subpath_t& subpath = mp_aln.subpath(i);
        auto& subpath_search_results = search_results[i];
        subpath_search_results.resize(subpath.path().mapping_size());
        for (size_t j = 0; j < subpath.path().mapping_size(); ++j) {
            // get the positions on paths that this mapping touches
            pos_t mapping_pos = make_pos_t(subpath.path().mapping(j).position());
            subpath_search_results[j] = nearest_offsets_in_paths(&graph, mapping_pos, 0, path_filter);
            // make sure that offsets are stored in increasing order
            for (pair<const path_handle_t, vector<pair<size_t, bool>>>& search_record : subpath_search_results[j]) {
                sort(search_record.second.begin(), search_record.second.end());
            }
#ifdef debug_mpaln_offsets
            cerr << "subpath " << i << ", mapping " << j << " path locations" << endl;
            for (const auto& pps : subpath_search_results[j]) {
                cerr << graph.get_path_name(pps.first) << endl;
                for (const auto& pp : pps.second) {
                    cerr << "\t" << pp.first << " " << pp.second << endl;
                }
            }
#endif
        }
    }
    
    path_positions_t return_val;
    
    // to keep track of whether we've already chosen a position on each path
    // earlier in the multipath alignment in either the forward or reverse pass
    vector<set<path_handle_t>> covered_fwd(mp_aln.subpath_size());
    vector<set<path_handle_t>> covered_rev(mp_aln.subpath_size());
    
    // forward pass looking for positions on the forward strand of paths
    for (size_t i = 0; i < mp_aln.subpath_size(); ++i) {
        const auto& subpath_search_results = search_results[i];
        for (size_t j = 0; j < subpath_search_results.size(); ++j) {
            for (const auto& path_pos : subpath_search_results[j]) {
                if (!covered_fwd[i].count(path_pos.first)) {
                    // we haven't already covered this path at an earlier position on the alignment
                    for (const auto& path_offset : path_pos.second) {
                        if (!path_offset.second) {
                            // there's a position on the forward strand of this path
                            return_val[path_pos.first].emplace_back(path_offset);
                            
                            // we're now covering this path for future search results
                            covered_fwd[i].insert(path_pos.first);
                            
#ifdef debug_mpaln_offsets
                            cerr << "found fwd pass pos, subpath " << i << ", mapping " << j << ", path " << graph.get_path_name(path_pos.first) << ", pos " << path_offset.first << " " << path_offset.second << endl;
#endif
                            
                            break;
                        }
                    }
                }
            }
        }
        
        // the following subpaths will be covered for any path that this
        // one is covered for
        for (auto n : mp_aln.subpath(i).next()) {
            auto& next_coverings = covered_fwd[n];
            for (auto path_handle : covered_fwd[i]) {
                next_coverings.insert(path_handle);
            }
        }
        for (const auto& c : mp_aln.subpath(i).connection()) {
            auto& next_coverings = covered_fwd[c.next()];
            for (auto path_handle : covered_fwd[i]) {
                next_coverings.insert(path_handle);
            }
        }
    }
    
    // now do a backward pass for the reverse strand of paths
    for (int64_t i = mp_aln.subpath_size() - 1; i >= 0; --i) {
        // find which paths are already covered in the reverse
        for (auto n : mp_aln.subpath(i).next()) {
            for (auto path_handle : covered_rev[n]) {
                covered_rev[i].insert(path_handle);
            }
        }
        for (const auto& c : mp_aln.subpath(i).connection()) {
            for (auto path_handle : covered_rev[c.next()]) {
                covered_rev[i].insert(path_handle);
            }
        }
        
        const auto& subpath_search_results = search_results[i];
        for (int64_t j = subpath_search_results.size() - 1; j >= 0; --j) {
            for (const auto& path_pos : subpath_search_results[j]) {
                if (!covered_rev[i].count(path_pos.first)) {
                    // we haven't already covered this path at an earlier position on the alignment
                    for (const auto& path_offset : path_pos.second) {
                        if (path_offset.second) {
                            // there's a position on the reverse strand of this path
                            auto mapping_len = mapping_from_length(mp_aln.subpath(i).path().mapping(j));
                            return_val[path_pos.first].emplace_back(path_offset.first - mapping_len,
                                                                    path_offset.second);
                            
#ifdef debug_mpaln_offsets
                            cerr << "found rev pass pos, subpath " << i << ", mapping " << j << ", path " << graph.get_path_name(path_pos.first) << ", pos " << path_offset.first - mapping_len << " " << path_offset.second << endl;
#endif
                            // we're now covering this path for future search results
                            covered_rev[i].insert(path_pos.first);
                            
                            break;
                        }
                    }
                }
            }
        }
    }
    
    return return_val;
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, int64_t search_limit, const std::function<bool(const path_handle_t&)>* path_filter) {
    annotate_with_path_positions(graph, aln, true, search_limit, path_filter);
}

void annotate_with_node_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, int64_t search_limit, const std::function<bool(const path_handle_t&)>* path_filter) {
    annotate_with_path_positions(graph, aln, false, search_limit, path_filter);
}

void annotate_with_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, bool just_min, int64_t search_limit, const std::function<bool(const path_handle_t&)>* path_filter) {
    if (!aln.refpos_size()) {
        // Get requested path positions
        unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions = alignment_path_offsets(graph, aln, just_min, false, search_limit, path_filter);
        // emit them in order of the path handle
        vector<path_handle_t> ordered;
        for (auto& path : positions) { ordered.push_back(path.first); }
        std::sort(ordered.begin(), ordered.end(), [](const path_handle_t& a, const path_handle_t& b) { return as_integer(a) < as_integer(b); });
        for (auto& path : ordered) {
            for (auto& p : positions[path]) {
                // Add each determined refpos
                Position* refpos = aln.add_refpos();
                refpos->set_name(graph.get_path_name(path));
                refpos->set_offset(p.first);
                refpos->set_is_reverse(p.second);
            }
        }
    }
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, vector<Alignment>& alns, int64_t search_limit, const std::function<bool(const path_handle_t&)>* path_filter) {
    for (auto& aln : alns) annotate_with_initial_path_positions(graph, aln, search_limit, path_filter);
}

}
}
