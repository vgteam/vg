#include "subgraph.hpp"

namespace vg {
namespace algorithms {

void expand_subgraph_by_steps(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t steps) {
    for (uint64_t i = 0; i < steps; ++i) {
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
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(x, h)) {
                        subgraph.create_edge(x, h);
                    }
                });
        }
    }
}

void expand_subgraph_to_node_count(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t node_count) {
    while (subgraph.get_node_count() < node_count && subgraph.get_node_count()) {
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
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (!subgraph.has_edge(x, h)) {
                        subgraph.create_edge(x, h);
                    }
                });
        }
    }
}

void expand_subgraph_by_length(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t length) {
    uint64_t accumulated_length = 0;
    while (accumulated_length < length) {
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
                        accumulated_length += subgraph.get_length(x);
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
}

void expand_subgraph_to_length(const HandleGraph& source, MutableHandleGraph& subgraph, uint64_t length) {
    uint64_t total_length = 0;
    subgraph.for_each_handle([&](const handle_t& h) {
            total_length += subgraph.get_length(h);
        });
    while (total_length < length) {
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
                        total_length += subgraph.get_length(x);
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

}
}
