#include "phase_unfolder.hpp"

#include <iostream>
#include <map>

namespace vg {

PhaseUnfolder::PhaseUnfolder(const xg::XG& xg_index, const gbwt::GBWT& gbwt_index, vg::id_t next_node) :
    xg_index(xg_index), gbwt_index(gbwt_index), next_node(next_node) {
}

void PhaseUnfolder::unfold(VG& graph, bool show_progress) {
    std::list<VG> components = this->complement_components(graph, show_progress);

    size_t haplotype_paths = 0;
    VG unfolded;
    for (VG& component : components) {
        haplotype_paths += this->unfold_component(component, graph, unfolded);
    }
    if (show_progress) {
        std::cerr << "Unfolded graph: "
                  << unfolded.node_count() << " nodes, " << unfolded.edge_count() << " edges on "
                  << haplotype_paths << " paths" << std::endl;
    }

    graph.extend(unfolded);
}

std::list<VG> PhaseUnfolder::complement_components(VG& graph, bool show_progress) {
    VG complement;

    // Add missing edges supported by XG paths.
    for (size_t path_rank = 1; path_rank <= this->xg_index.max_path_rank(); path_rank++) {
        const xg::XGPath& path = this->xg_index.get_path(this->xg_index.path_name(path_rank));
        size_t path_length = path.ids.size();
        if (path_length == 0) {
            continue;
        }
        gbwt::node_type prev = gbwt::Node::encode(path.node(0), path.is_reverse(0));
        for (size_t i = 1; i < path_length; i++) {
            gbwt::node_type curr = gbwt::Node::encode(path.node(i), path.is_reverse(i));
            Edge candidate = xg::make_edge(gbwt::Node::id(prev), gbwt::Node::is_reverse(prev),
                                           gbwt::Node::id(curr), gbwt::Node::is_reverse(curr));
            if (!graph.has_edge(candidate)) {
                complement.add_node(this->xg_index.node(candidate.from()));
                complement.add_node(this->xg_index.node(candidate.to()));
                complement.add_edge(candidate);
            }
            prev = curr;
        }
    }

    // Add missing edges supported by GBWT threads.
    for (gbwt::comp_type comp = 1; comp < this->gbwt_index.effective(); comp++) {
        gbwt::node_type gbwt_node = this->gbwt_index.toNode(comp);
        std::vector<gbwt::edge_type> outgoing = this->gbwt_index.edges(gbwt_node);
        for (gbwt::edge_type outedge : outgoing) {
            if (outedge.first == gbwt::ENDMARKER) {
                continue;
            }
            Edge candidate = xg::make_edge(gbwt::Node::id(gbwt_node), gbwt::Node::is_reverse(gbwt_node),
                                           gbwt::Node::id(outedge.first), gbwt::Node::is_reverse(outedge.first));
            if (!graph.has_edge(candidate)) {
                complement.add_node(this->xg_index.node(candidate.from()));
                complement.add_node(this->xg_index.node(candidate.to()));
                complement.add_edge(candidate);
            }
        }
    }

    std::list<VG> components;
    complement.disjoint_subgraphs(components);
    if (show_progress) {
        std::cerr << "Complement graph: "
                  << complement.node_count() << " nodes, " << complement.edge_count() << " edges in "
                  << components.size() << " components" << std::endl;
    }
    return components;
}

size_t PhaseUnfolder::unfold_component(VG& component, VG& graph, VG& unfolded) {
    // Find the border nodes shared between the component and the graph.
    component.for_each_node([&](Node* node) {
       if (graph.has_node(node->id())) {
           this->border.insert(node->id());
       }
    });

    // Generate the paths starting from each border node.
    for (vg::id_t start_node : this->border) {
        this->generate_paths(component, start_node);
        this->generate_threads(component, start_node);
    }
    size_t haplotype_paths = this->paths.size();

    // Unfold the generated paths. We merge duplicated nodes by shared
    // prefixes in the first half of the path and by shared suffixes in
    // the second half. Needless duplication would otherwise make GCSA2
    // index construction too expensive.
    // TODO: We need to store the mapping.
    std::map<path_type, vg::id_t> node_by_prefix, node_by_suffix;
    for (const path_type& path : this->paths) {
        Node prev = this->xg_index.node(gbwt::Node::id(path.front()));
        unfolded.add_node(prev);
        for (size_t i = 1; i < path.size(); i++) {
            Node curr = this->xg_index.node(gbwt::Node::id(path[i]));
            bool is_prefix = (i < (path.size() + 1) / 2);
            bool is_suffix = !is_prefix & (i + 1 < path.size());
            if (is_prefix) {
                auto iter = node_by_prefix.emplace(path_type(path.begin(), path.begin() + i + 1), 0).first;
                if (iter->second == 0) {      // No cached node with the same prefix.
                    iter->second = this->get_id();
                }
                curr.set_id(iter->second);
            } else if (is_suffix) {
                auto iter = node_by_suffix.emplace(path_type(path.begin() + i, path.end()), 0).first;
                if (iter->second == 0) {      // No cached node with the same suffix.
                    iter->second = this->get_id();
                }
                curr.set_id(iter->second);
            }
            unfolded.add_node(curr);
            Edge edge = xg::make_edge(prev.id(), gbwt::Node::is_reverse(path[i - 1]),
                                      curr.id(), gbwt::Node::is_reverse(path[i]));
            unfolded.add_edge(edge);
            prev = curr;
        }
    }

    this->border.clear();
    this->paths.clear();
    return haplotype_paths;
}

void PhaseUnfolder::generate_paths(VG& component, vg::id_t from) {
    static int component_id = 0;
    for (size_t path_rank = 1; path_rank <= this->xg_index.max_path_rank(); path_rank++) {
        const xg::XGPath& path = this->xg_index.get_path(this->xg_index.path_name(path_rank));
        size_t path_length = path.ids.size();

        std::vector<size_t> occurrences = this->xg_index.node_ranks_in_path(from, path_rank);
        for (size_t occurrence : occurrences) {
            gbwt::node_type prev = gbwt::Node::encode(path.node(occurrence), path.is_reverse(occurrence));
            path_type buffer { prev };
            for (size_t i = occurrence + 1; i < path_length; i++) {
                gbwt::node_type curr = gbwt::Node::encode(path.node(i), path.is_reverse(i));
                Edge candidate = xg::make_edge(gbwt::Node::id(prev), gbwt::Node::is_reverse(prev),
                                               gbwt::Node::id(curr), gbwt::Node::is_reverse(curr));
                if (!component.has_edge(candidate)) {
                    break;
                }
                buffer.push_back(curr);
                if (this->border.find(gbwt::Node::id(curr)) != this->border.end()) {
                    this->insert_path(buffer);
                }
                prev = curr;
            }
        }
    }
}

void PhaseUnfolder::generate_threads(VG& component, vg::id_t from) {
    this->create_state(from, false);
    this->create_state(from, true);

    while (!this->states.empty()) {
        state_type state = this->states.top(); this->states.pop();
        vg::id_t node = gbwt::Node::id(state.second.back());
        bool is_reverse = gbwt::Node::is_reverse(state.second.back());

        std::vector<Edge*> edges = component.edges_of(component.get_node(node));
        for (Edge* edge : edges) {
            if (edge->from() == node && edge->from_start() == is_reverse) {
                this->extend_state(state, edge->to(), edge->to_end());
            }
            else if (edge->to() == node && edge->to_end() != is_reverse) {
                this->extend_state(state, edge->from(), !edge->from_start());
            }
        }

        if (this->border.find(node) != this->border.end()) {
            this->insert_path(state.second);
        }
    }
}

void PhaseUnfolder::create_state(vg::id_t node, bool is_reverse) {
    search_type search = this->gbwt_index.find(gbwt::Node::encode(node, is_reverse));
    if (search.empty()) {
        return;
    }
    this->states.push(std::make_pair(search, path_type {search.node}));
}

void PhaseUnfolder::extend_state(state_type state, vg::id_t node, bool is_reverse) {
    state.first = this->gbwt_index.extend(state.first, gbwt::Node::encode(node, is_reverse));
    if (state.first.empty()) {
        return;
    }
    state.second.push_back(state.first.node);
    this->states.push(state);
}

void PhaseUnfolder::insert_path(const path_type& path) {
    if (path.size() < 2) {
        return;
    }
    path_type reverse_complement(path.size(), 0);
    for (size_t i = 0; i < path.size(); i++) {
        reverse_complement[path.size() - 1 - i] = gbwt::Node::reverse(path[i]);
    }
    this->paths.insert(std::min(path, reverse_complement));
}

vg::id_t PhaseUnfolder::get_id() {
    return this->next_node++;
}

} 
