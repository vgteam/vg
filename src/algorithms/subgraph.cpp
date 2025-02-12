#include "subgraph.hpp"
#include "../path.hpp"
#include "../crash.hpp"

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
        if (!subgraph.has_node(i) && source.has_node(i)) {
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

void add_subpaths_to_subgraph(const PathPositionHandleGraph& source, MutablePathHandleGraph& subgraph) {

    // We want to organize all visits by base path. This key type holds the
    // sense, sample and locus names, haplotype, and phase block.
    using base_metadata_t = std::tuple<PathSense, string, string, size_t, size_t>;

    // This stores, for each source graph base path, for each start offset, the handle at that offset on the path.
    std::unordered_map<base_metadata_t, std::map<uint64_t, handle_t> > subpaths;

    // This stores information about base paths that don't have subranges, and
    // their full lengths and circularity flags, so we can avoid generating new
    // subrange metadata when we just have all of a path.
    std::unordered_map<base_metadata_t, std::pair<size_t, bool>> full_path_info;

    subgraph.for_each_handle([&](const handle_t& h) {
            handlegraph::nid_t id = subgraph.get_id(h);
            if (source.has_node(id)) {
                handle_t handle = source.get_handle(id);
                source.for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                        path_handle_t path = source.get_path_handle_of_step(step);
                        // Figure out the base path this visit is on
                        base_metadata_t key = {source.get_sense(path), source.get_sample_name(path), source.get_locus_name(path), source.get_haplotype(path), source.get_phase_block(path)};
                        // Figure out the subrange of the base path it is relative to
                        subrange_t path_subrange = source.get_subrange(path);
                        uint64_t visit_offset = pos;
                        if (path_subrange != PathMetadata::NO_SUBRANGE) {
                            // If we have the position relative to a subrange, adjust by that subrange's offset.
                            visit_offset += path_subrange.first;
                        }
                        subpaths[key][visit_offset] = is_rev ? subgraph.flip(h) : h;
                        
                        if (path_subrange == PathMetadata::NO_SUBRANGE) {
                            // There's no subrange set, so this path is full-length in the source graph.
                            // See if we know of this path as a full-length path or not
                            auto it = full_path_info.find(key);
                            if (it == full_path_info.end()) {
                                // We haven't recorded its length and circularity yet, so do it.
                                full_path_info.emplace_hint(it, key, std::make_pair(source.get_path_length(path), source.get_is_circular(path)));
                            }
                        }
                        return true;
                    });
            }
        });

    for (auto& base_and_visits : subpaths) {
        // For each base path
        const base_metadata_t& base_path_metadata = base_and_visits.first;
        const auto& start_to_handle = base_and_visits.second;
        // If we didn't put anything in the visit collection, it shouldn't be here.
        crash_unless(!start_to_handle.empty());

        // We're going to walk over all the visits and find contiguous runs
        auto run_start = start_to_handle.begin();
        auto run_end = run_start;
        size_t start_coordinate = run_start->first;
        while (run_end != start_to_handle.end()) {
            // Until we run out of runs
            // Figure out where this node ends on the path
            size_t stop_coordinate = run_end->first + subgraph.get_length(run_end->second);
            
            // Look ahead
            ++run_end;

            if (run_end != start_to_handle.end() && run_end->first == stop_coordinate) {
                // The next visit is still contiguous, so advance.
                continue;
            }

            // Otherwise we've reached a break in continuity. We have a
            // contiguous run from run_start to run_end, visiting the subrange
            // start_coordinate to stop_coordinate.

            // Find out if we cover a full source graph path.
            subrange_t run_subrange = {start_coordinate, stop_coordinate};
            bool is_circular = false;
            if (start_coordinate == 0) {
                // We might be a full path
                auto found_length_and_circularity = full_path_info.find(base_path_metadata);
                if (found_length_and_circularity != full_path_info.end() && found_length_and_circularity->second.first == stop_coordinate) {
                    // We are a full path
                    run_subrange = PathMetadata::NO_SUBRANGE;
                    // We can be circular.
                    is_circular = found_length_and_circularity->second.second;
                }
            }
            
            // Make a path with all the metadata
            path_handle_t new_path = subgraph.create_path(
                std::get<0>(base_path_metadata),
                std::get<1>(base_path_metadata),
                std::get<2>(base_path_metadata),
                std::get<3>(base_path_metadata),
                std::get<4>(base_path_metadata),
                run_subrange,
                is_circular
            );

            for (auto it = run_start; it != run_end; ++it) {
                // Copy the path's visits
                subgraph.append_step(new_path, it->second);
            }

            // Set up the next subpath.
            // Set where it starts.
            run_start = run_end;
            if (run_start != start_to_handle.end()) {
                // And if it will exist, set its start coordinate.
                start_coordinate = run_start->first;
            }
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
