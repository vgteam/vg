#include "subgraph.hpp"
#include "../path.hpp"

namespace vg {
namespace algorithms {

void expand_subgraph_by_steps(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& steps, bool forward_only) {
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
    for (uint64_t i = 0; i < steps && !curr_handles.empty(); ++i) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                        next_handles.push_back(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (source.get_is_reverse(c)) {
                        x = subgraph.flip(x);
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(x, h)) {
                            subgraph.create_edge(x, h);
                        }
                    });
            }
        }
        curr_handles = std::move(next_handles);
    }
    add_connecting_edges_to_subgraph(source, subgraph);
}

void expand_subgraph_to_node_count(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& node_count, bool forward_only) {
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
    while (subgraph.get_node_count() < node_count && subgraph.get_node_count()) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                        next_handles.push_back(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (source.get_is_reverse(c)) {
                        x = subgraph.flip(x);
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(x, h)) {
                            subgraph.create_edge(x, h);
                        }
                    });
            }
        }
        curr_handles = std::move(next_handles);
    }
    add_connecting_edges_to_subgraph(source, subgraph);
}

void expand_subgraph_by_length(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& length, bool forward_only) {
    uint64_t accumulated_length = 0;
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
    while (accumulated_length < length && !curr_handles.empty()) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                        next_handles.push_back(x);
                        accumulated_length += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (source.get_is_reverse(c)) {
                        x = subgraph.flip(x);
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                            next_handles.push_back(x);
                            accumulated_length += subgraph.get_length(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(x, h)) {
                            subgraph.create_edge(x, h);
                        }
                    });
            }
        }
        curr_handles = std::move(next_handles);
    }
    add_connecting_edges_to_subgraph(source, subgraph);
}

void expand_subgraph_to_length(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& length, bool forward_only) {
    uint64_t total_length = 0;
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            total_length += subgraph.get_length(h);
            curr_handles.push_back(h);
        });
    while (total_length < length && !curr_handles.empty()) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                        next_handles.push_back(x);
                        total_length += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (source.get_is_reverse(c)) {
                        x = subgraph.flip(x);
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                            next_handles.push_back(x);
                            total_length += subgraph.get_length(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(x, h)) {
                            subgraph.create_edge(x, h);
                        }
                    });
            }
        }
        curr_handles = std::move(next_handles);
    }
    add_connecting_edges_to_subgraph(source, subgraph);
}

/// expand the context around a single handle position
void extract_context(const HandleGraph& source, MutableHandleGraph& subgraph, const handle_t& handle, const uint64_t& offset, const uint64_t& length, bool fwd, bool rev) {
    uint64_t total_length_fwd = source.get_length(handle)-offset;
    uint64_t total_length_rev = offset;
    uint64_t get_fwd = fwd && !rev ? length : length/2;
    uint64_t get_rev = !fwd && rev ? length : length/2;
    if (!subgraph.has_node(source.get_id(handle))) {
        subgraph.create_handle(source.get_sequence(source.get_is_reverse(handle)?source.flip(handle):handle), source.get_id(handle));
    }
    bool extended = true;
    while (extended && (total_length_fwd < get_fwd || total_length_rev < get_rev)) {
        std::vector<handle_t> curr_handles;
        subgraph.for_each_handle([&](const handle_t& h) {
                curr_handles.push_back(h);
            });
        extended = false;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            if (total_length_fwd < get_fwd) {
                source.follow_edges(old_h, false, [&](const handle_t& c) {
                        if (total_length_fwd < get_fwd) {
                            handle_t x;
                            if (!subgraph.has_node(source.get_id(c))) {
                                x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                                total_length_fwd += subgraph.get_length(x);
                                extended = true;
                            } else {
                                x = subgraph.get_handle(source.get_id(c));
                            }
                            if (source.get_is_reverse(c)) {
                                x = subgraph.flip(x);
                            }
                            if (!subgraph.has_edge(h, x)) {
                                subgraph.create_edge(h, x);
                            }
                        }
                    });
            }
            if (total_length_rev < get_rev) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        if (total_length_rev < get_rev) {
                            handle_t x;
                            if (!subgraph.has_node(source.get_id(c))) {
                                x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                                total_length_rev += subgraph.get_length(x);
                                extended = true;
                            } else {
                                x = subgraph.get_handle(source.get_id(c));
                            }
                            if (source.get_is_reverse(c)) {
                                x = subgraph.flip(x);
                            }
                            if (!subgraph.has_edge(x, h)) {
                                subgraph.create_edge(x, h);
                            }
                        }
                    });
            }
        }
    }
    add_connecting_edges_to_subgraph(source, subgraph);
}

void extract_id_range(const HandleGraph& source, const nid_t& id1, const nid_t& id2, MutableHandleGraph& subgraph) {
    for (nid_t i = id1; i <= id2; ++i) {
        if (!subgraph.has_node(i)) {
            subgraph.create_handle(source.get_sequence(source.get_handle(i)), i);
        }
    }
}

void extract_path_range(const PathPositionHandleGraph& source, path_handle_t path_handle, int64_t start, int64_t end,
                        MutableHandleGraph& subgraph) {
    step_handle_t start_step = source.get_step_at_position(path_handle, start);
    size_t start_position = source.get_position_of_step(start_step);
    size_t size_needed = end < 0 ? numeric_limits<size_t>::max() : end - start + 1 + start - start_position;
    size_t running_length = 0;
    
    for (step_handle_t cur_step = start_step; cur_step != source.path_end(path_handle) && running_length < size_needed;
         cur_step = source.get_next_step(cur_step)) {
        handle_t cur_handle = source.get_handle_of_step(cur_step);
        subgraph.create_handle(source.get_sequence(cur_handle), source.get_id(cur_handle));
        if (cur_step != start_step) {
            handle_t prev_handle = source.get_handle_of_step(source.get_previous_step(cur_step));
            subgraph.create_edge(subgraph.get_handle(source.get_id(prev_handle), source.get_is_reverse(prev_handle)),
                                  subgraph.get_handle(source.get_id(cur_handle), source.get_is_reverse(cur_handle)));
        }
        running_length += source.get_length(cur_handle);
    }
}

/// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
/// based on their order in the path position index provided by the source graph
/// will clear any path found in both graphs before writing the new steps into it
/// if subpath_naming is true, a suffix will be added to each path in the subgraph denoting its offset
/// in the source graph (unless the subpath was not cut up at all)
void add_subpaths_to_subgraph(const PathPositionHandleGraph& source, MutablePathHandleGraph& subgraph,
                              bool subpath_naming) {
    std::unordered_map<std::string, std::map<uint64_t, handle_t> > subpaths;
    subgraph.for_each_handle([&](const handle_t& h) {
            handlegraph::nid_t id = subgraph.get_id(h);
            if (source.has_node(id)) {
                handle_t handle = source.get_handle(id);
                source.for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                        path_handle_t path = source.get_path_handle_of_step(step);
                        std::string path_name = source.get_path_name(path);
                        subpaths[path_name][pos] = is_rev ? subgraph.flip(h) : h;
                        return true;
                    });
            }
        });

    function<path_handle_t(const path_handle_t&, bool, size_t)> new_subpath =
        [&source, &subgraph](const path_handle_t& source_path_handle, bool is_circular, size_t subpath_offset) {
        
        subrange_t subrange = source.get_subrange(source_path_handle);

        // TODO: Are we really supposed to ignore the original path's start offset?
        if (subrange == PathMetadata::NO_SUBRANGE) {
            subrange.first = subpath_offset;
        } else {
            subrange.first += subpath_offset;
        }
        // TODO: Are we really supposed to ignore the original path's start
        // offset? We used to clobber it completely with subpath_offset.
        subrange.second = PathMetadata::NO_END_POSITION;

        PathSense sense = source.get_sense(source_path_handle);
        std::string sample = source.get_sample_name(source_path_handle);
        std::string locus = source.get_locus_name(source_path_handle);
        size_t haplotype = source.get_haplotype(source_path_handle);

        // TODO: We're supposed to clobber existing paths, so there might be an
        // existing path like the one we're trying to make here. But we don't
        // have a has_path or efficient get_path by full metadata. So we make
        // the name twice.
        string subpath_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, subrange);
        if (subgraph.has_path(subpath_name)) {
            subgraph.destroy_path(subgraph.get_path_handle(subpath_name));
        }
        return subgraph.create_path(sense, sample, locus, haplotype, subrange, is_circular);
    };

    for (auto& subpath : subpaths) {
        const std::string& path_name = subpath.first;
        path_handle_t source_path_handle = source.get_path_handle(path_name);
        // destroy the path if it exists
        if (subgraph.has_path(path_name)) {
            subgraph.destroy_path(subgraph.get_path_handle(path_name));
        }
        // create a new path.  give it a subpath name if the flag's on and its smaller than original
        path_handle_t path;
        if (!subpath_naming || subpath.second.size() == source.get_step_count(source_path_handle) ||
            subpath.second.empty()) {
            path = subgraph.create_path_handle(path_name, source.get_is_circular(source_path_handle));
        } else {
            path = new_subpath(source_path_handle, source.get_is_circular(source_path_handle), subpath.second.begin()->first);
        }
        for (auto p = subpath.second.begin(); p != subpath.second.end(); ++p) {
            const handle_t& handle = p->second;
            if (p != subpath.second.begin() && subpath_naming) {
                auto prev = p;
                --prev;
                const handle_t& prev_handle = prev->second;
                // distance from map
                size_t delta = p->first - prev->first;
                // what the distance should be if they're contiguous depends on relative orienations
                size_t cont_delta = subgraph.get_length(prev_handle);
                if (delta != cont_delta) {
                    // we have a discontinuity!  we'll make a new path can continue from there
                    assert(subgraph.get_step_count(path) > 0);
                    path = new_subpath(source_path_handle, subgraph.get_is_circular(path), p->first);
                }
            }
            //fill in the path information
            subgraph.append_step(path, handle);
        }
    }
}

/// We can accumulate a subgraph without accumulating all the edges between its nodes
/// this helper ensures that we get the full set
void add_connecting_edges_to_subgraph(const HandleGraph& source, MutableHandleGraph& subgraph) {
    subgraph.for_each_handle([&](const handle_t& handle) {
            nid_t id = subgraph.get_id(handle);
            handle_t source_handle = source.get_handle(id, subgraph.get_is_reverse(handle));
            source.follow_edges(source_handle, false, [&](const handle_t& next) {
                    nid_t next_id = source.get_id(next);
                    if (subgraph.has_node(next_id)) {
                        handle_t subgraph_next = subgraph.get_handle(next_id, source.get_is_reverse(next));
                        if (!subgraph.has_edge(handle, subgraph_next)) {
                            subgraph.create_edge(handle, subgraph_next);
                        }
                    }
                });
            source.follow_edges(source_handle, true, [&](const handle_t& prev) {
                    nid_t prev_id = source.get_id(prev);
                    if (subgraph.has_node(prev_id)) {
                        handle_t subgraph_prev = subgraph.get_handle(prev_id, source.get_is_reverse(prev));
                        if (!subgraph.has_edge(subgraph_prev, handle)) {
                            subgraph.create_edge(subgraph_prev, handle);
                        }
                    }
                });
        });
}

}
}
