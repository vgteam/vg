#include "phase_unfolder.hpp"

#include <iostream>

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
    for (gbwt::comp_type comp = 1; comp < this->gbwt_index.effective(); comp++) {
        gbwt::node_type gbwt_node = this->gbwt_index.toNode(comp);
        std::vector<gbwt::edge_type> outgoing = this->gbwt_index.edges(gbwt_node);
        for (gbwt::edge_type outedge : outgoing) {
            if (outedge.first == gbwt::ENDMARKER) {
                continue;
            }
            Edge candidate;
            candidate.set_from(gbwt::Node::id(gbwt_node));
            candidate.set_to(gbwt::Node::id(outedge.first));
            candidate.set_from_start(gbwt::Node::is_reverse(gbwt_node));
            candidate.set_to_end(gbwt::Node::is_reverse(outedge.first));
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
    }
    size_t haplotype_paths = this->paths.size();

    // Unfold the generated paths.
    // TODO: We should merge shared prefixes and/or suffixes.
    for (const path_type& path : this->paths) {
        Node prev = this->xg_index.node(gbwt::Node::id(path.front()));
        unfolded.add_node(prev);
        for (size_t i = 1; i < path.size(); i++) {
            Node node = this->xg_index.node(gbwt::Node::id(path[i]));
            if (i + 1 < path.size()) {
                node.set_id(this->next_node); this->next_node++; // TODO: We need to store the mapping.
            }
            unfolded.add_node(node);
            Edge edge;
            edge.set_from(prev.id());
            edge.set_to(node.id());
            edge.set_from_start(gbwt::Node::is_reverse(path[i - 1]));
            edge.set_to_end(gbwt::Node::is_reverse(path[i]));
            unfolded.add_edge(edge);
            prev = node;
        }
    }

    this->border.clear();
    this->paths.clear();
    return haplotype_paths;
}

void PhaseUnfolder::generate_paths(VG& component, vg::id_t from) {
    this->create_state(from, false);
    this->create_state(from, true);

    while (!this->states.empty()) {
        state_type state = this->states.top(); this->states.pop();
        vg::id_t node = gbwt::Node::id(state.second.back());
        bool is_reverse = gbwt::Node::is_reverse(state.second.back());

        std::vector<Edge*> edges = component.edges_of(component.get_node(node));
        bool has_successors = false;
        for (Edge* edge : edges) {
            if (edge->from() == node && edge->from_start() == is_reverse) {
                has_successors |= this->extend_state(state, edge->to(), edge->to_end());
            }
            else if (edge->to() == node && edge->to_end() != is_reverse) {
                has_successors |= this->extend_state(state, edge->from(), !edge->from_start());
            }
        }

        if (!has_successors || this->border.find(node) != this->border.end()) {
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

bool PhaseUnfolder::extend_state(state_type state, vg::id_t node, bool is_reverse) {
    state.first = this->gbwt_index.extend(state.first, gbwt::Node::encode(node, is_reverse));
    if (state.first.empty()) {
        return false;
    }
    state.second.push_back(state.first.node);
    this->states.push(state);
    return true;
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

} 
