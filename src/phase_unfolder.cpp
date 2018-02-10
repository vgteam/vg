#include "phase_unfolder.hpp"

#include <iostream>
#include <map>

namespace vg {

PhaseUnfolder::PhaseUnfolder(const xg::XG& xg_index, const gbwt::GBWT& gbwt_index, vg::id_t next_node) :
    xg_index(xg_index), gbwt_index(gbwt_index), mapping(next_node) {
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

void PhaseUnfolder::restore_paths(VG& graph, bool show_progress) {

    for (size_t path_rank = 1; path_rank <= this->xg_index.max_path_rank(); path_rank++) {
        const xg::XGPath& path = this->xg_index.get_path(this->xg_index.path_name(path_rank));
        size_t path_length = path.ids.size();
        if (path_length == 0) {
            continue;
        }

        gbwt::node_type prev = gbwt::Node::encode(path.node(0), path.is_reverse(0));
        for (size_t i = 1; i < path_length; i++) {
            gbwt::node_type curr = gbwt::Node::encode(path.node(i), path.is_reverse(i));
            Edge candidate = make_edge(prev, curr);
            if (!graph.has_edge(candidate)) {
                graph.add_node(this->xg_index.node(candidate.from()));
                graph.add_node(this->xg_index.node(candidate.to()));
                graph.add_edge(candidate);
            }
            prev = curr;
        }
    }

    if (show_progress) {
        std::cerr << "Restored graph: " << graph.node_count() << " nodes, " << graph.edge_count() << " edges" << std::endl;
    }
}

void PhaseUnfolder::write_mapping(const std::string& filename) const {
    std::ofstream out(filename, std::ios_base::binary);
    if (!out) {
        std::cerr << "[PhaseUnfolder]: cannot create mapping file " << filename << std::endl;
        return;
    }
    this->mapping.serialize(out);
    out.close();
}

void PhaseUnfolder::read_mapping(const std::string& filename) {
    std::ifstream in(filename, std::ios_base::binary);
    if (!in) {
        std::cerr << "[PhaseUnfolder]: cannot open mapping file " << filename << std::endl;
        return;
    }
    this->mapping.load(in);
    in.close();
}

vg::id_t PhaseUnfolder::get_mapping(vg::id_t node) const {
    return this->mapping(node);
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
            Edge candidate = make_edge(prev, curr);
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
            Edge candidate = make_edge(gbwt_node, outedge.first);
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

    auto insert_node = [&](gbwt::node_type node) {
        Node temp = this->xg_index.node(this->get_mapping(gbwt::Node::id(node)));
        temp.set_id(gbwt::Node::id(node));
        unfolded.add_node(temp);
    };

    // Create the unfolded component from the tries.
    for (auto mapping : this->prefixes) {
        gbwt::node_type from = mapping.first.first, to = mapping.second;
        insert_node(from);
        insert_node(to);
        unfolded.add_edge(make_edge(from, to));
    }
    for (auto mapping : this->suffixes) {
        gbwt::node_type from = mapping.second, to = mapping.first.second;
        insert_node(from);
        insert_node(to);
        unfolded.add_edge(make_edge(from, to));
    }
    for (auto edge : this->crossing_edges) {
        insert_node(edge.first);
        insert_node(edge.second);
        unfolded.add_edge(make_edge(edge.first, edge.second));
    }

    size_t haplotype_paths = this->crossing_edges.size();
    this->border.clear();
    this->prefixes.clear();
    this->suffixes.clear();
    this->crossing_edges.clear();
    return haplotype_paths;
}

void PhaseUnfolder::generate_paths(VG& component, vg::id_t from) {

    for (size_t path_rank = 1; path_rank <= this->xg_index.max_path_rank(); path_rank++) {
        const xg::XGPath& path = this->xg_index.get_path(this->xg_index.path_name(path_rank));
        size_t path_length = path.ids.size();

        std::vector<size_t> occurrences = this->xg_index.node_ranks_in_path(from, path_rank);
        for (size_t occurrence : occurrences) {
            // Forward.
            {
                gbwt::node_type prev = gbwt::Node::encode(path.node(occurrence), path.is_reverse(occurrence));
                path_type buffer { prev };
                for (size_t i = occurrence + 1; i < path_length; i++) {
                    gbwt::node_type curr = gbwt::Node::encode(path.node(i), path.is_reverse(i));
                    Edge candidate = make_edge(prev, curr);
                    if (!component.has_edge(candidate)) {
                        break;  // Found a maximal path.
                    }
                    buffer.push_back(curr);
                    if (this->border.find(gbwt::Node::id(curr)) != this->border.end()) {
                        break;  // Found a border-to-border path.
                    }
                    prev = curr;
                }
                this->insert_path(buffer);
            }

            // Backward.
            {
                gbwt::node_type prev = gbwt::Node::encode(path.node(occurrence), !path.is_reverse(occurrence));
                path_type buffer { prev };
                bool found_border = false;
                for (size_t i = occurrence; i > 0 ; i--) {
                    gbwt::node_type curr = gbwt::Node::encode(path.node(i - 1), !path.is_reverse(i - 1));
                    Edge candidate = make_edge(prev, curr);
                    if (!component.has_edge(candidate)) {
                        break;  // Found a maximal path.
                    }
                    buffer.push_back(curr);
                    if (this->border.find(gbwt::Node::id(curr)) != this->border.end()) {
                        break;  // Found a border-to-border path.
                    }
                    prev = curr;
                }
                this->insert_path(buffer);
            }
        }
    }
}

void PhaseUnfolder::generate_threads(VG& component, vg::id_t from) {
    this->create_state(from, false);
    this->create_state(from, true);

    while (!this->states.empty()) {
        state_type state = this->states.top(); this->states.pop();
        vg::id_t node = gbwt::Node::id(state.first.node);
        bool is_reverse = gbwt::Node::is_reverse(state.first.node);

        if (state.second.size() >= 2 && this->border.find(node) != this->border.end()) {
            this->insert_path(state.second);    // Border-to-border path.
            continue;
        }

        std::vector<Edge*> edges = component.edges_of(component.get_node(node));
        bool was_extended = false;
        for (Edge* edge : edges) {
            if (edge->from() == node && edge->from_start() == is_reverse) {
                was_extended |= this->extend_state(state, edge->to(), edge->to_end());
            }
            else if (edge->to() == node && edge->to_end() != is_reverse) {
                was_extended |= this->extend_state(state, edge->from(), !edge->from_start());
            }
        }

        if (!was_extended) {
            this->insert_path(state.second);    // Maximal path.
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
    const path_type& to_insert = std::min(path, reverse_complement);

    /*
      Break the path in half. For each prefix / suffix, check if we already
      have a mapping for the next node. If not, create a new duplicate of
      the node and insert the mapping into the corresponding trie. Finally
      create a crossing edge between the full prefix and the full suffix.
    */

    // Prefixes.
    gbwt::node_type from = to_insert.front();
    for (size_t i = 1; i < (to_insert.size() + 1) / 2; i++) {
        auto iter = this->prefixes.emplace(std::make_pair(from, to_insert[i]), 0).first;
        if (iter->second == 0) {
            gbwt::size_type new_id = this->mapping.insert(gbwt::Node::id(to_insert[i]));
            iter->second = gbwt::Node::encode(new_id, gbwt::Node::is_reverse(to_insert[i]));
        }
        from = iter->second;
    }

    // Suffixes.
    gbwt::node_type to = to_insert.back();
    for (size_t i = to_insert.size() - 2; i >= (to_insert.size() + 1) / 2; i--) {
        auto iter = this->suffixes.emplace(std::make_pair(to_insert[i], to), 0).first;
        if (iter->second == 0) {
            gbwt::size_type new_id = this->mapping.insert(gbwt::Node::id(to_insert[i]));
            iter->second = gbwt::Node::encode(new_id, gbwt::Node::is_reverse(to_insert[i]));
        }
        to = iter->second;
    }

    // Crossing edge.
    this->crossing_edges.insert(std::make_pair(from, to));
}

} 
