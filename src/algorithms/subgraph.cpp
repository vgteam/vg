#include "subgraph.hpp"

namespace vg {
namespace algorithms {

void expand_subgraph_by_steps(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t steps, bool forward_only) {
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
    for (uint64_t i = 0; i < steps; ++i) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        next_handles.push_back(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
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

void expand_subgraph_to_node_count(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t node_count, bool forward_only) {
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
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        next_handles.push_back(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
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

void expand_subgraph_by_length(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t length, bool forward_only) {
    uint64_t accumulated_length = 0;
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
    while (accumulated_length < length) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        next_handles.push_back(x);
                        accumulated_length += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                            next_handles.push_back(x);
                            accumulated_length += subgraph.get_length(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
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

void expand_subgraph_to_length(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t length, bool forward_only) {
    uint64_t total_length = 0;
    std::vector<handle_t> curr_handles;
    subgraph.for_each_handle([&](const handle_t& h) {
            total_length += subgraph.get_length(h);
            curr_handles.push_back(h);
        });
    while (total_length < length) {
        std::vector<handle_t> next_handles;
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        next_handles.push_back(x);
                        total_length += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            if (!forward_only) {
                source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                            next_handles.push_back(x);
                            total_length += subgraph.get_length(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
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
void extract_context(const HandleGraph& source, MutableHandleGraph& subgraph, const handle_t& handle, const uint64_t& offset, const uint64_t& length) {
    uint64_t total_length_fwd = source.get_length(handle)-offset;
    uint64_t total_length_rev = offset;
    if (!subgraph.has_node(source.get_id(handle))) {
        subgraph.create_handle(source.get_sequence(handle), source.get_id(handle));
    }
    while (total_length_fwd < length/2 || total_length_rev < length/2) {
        std::vector<handle_t> curr_handles;
        subgraph.for_each_handle([&](const handle_t& h) {
                curr_handles.push_back(h);
            });
        for (auto& h : curr_handles) {
            handle_t old_h = source.get_handle(subgraph.get_id(h));
            source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        total_length_fwd += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
            source.follow_edges(old_h, true, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(c), source.get_id(c));
                        total_length_rev += subgraph.get_length(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(x, h)) {
                        subgraph.create_edge(x, h);
                    }
                });
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

/// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
/// based on their order in the path position index provided by the source graph
/// will clear any path found in both graphs before writing the new steps into it
void add_subpaths_to_subgraph(const PathPositionHandleGraph& source, MutablePathHandleGraph& subgraph) {
    std::unordered_map<std::string, std::map<uint64_t, handle_t> > subpaths;
    subgraph.for_each_handle([&](const handle_t& h) {
            handlegraph::nid_t id = subgraph.get_id(h);
            if (source.has_node(id)) {
                handle_t handle = source.get_handle(id);
                source.for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                        path_handle_t path = source.get_path_handle_of_step(step);
                        std::string path_name = source.get_path_name(path);
                        subpaths[path_name][pos] = h;
                        return true;
                    });
            }
        });
    for (auto& subpath : subpaths) {
        const std::string& path_name = subpath.first;
        // destroy the path if it exists
        if (subgraph.has_path(path_name)) {
            subgraph.destroy_path(subgraph.get_path_handle(path_name));
        }
        // fill in the path information
        path_handle_t path = subgraph.create_path_handle(path_name);
        for (auto& p : subpath.second) {
            const handle_t& handle = p.second;
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
