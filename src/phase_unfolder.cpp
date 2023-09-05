#include "phase_unfolder.hpp"
#include "progress_bar.hpp"
#include "algorithms/disjoint_components.hpp"

#include <cassert>
#include <iostream>
#include <map>
#include <set>

namespace vg {

PhaseUnfolder::PhaseUnfolder(const PathHandleGraph& path_graph, const gbwt::GBWT& gbwt_index, vg::id_t next_node) :
    path_graph(path_graph), gbwt_index(gbwt_index), mapping(next_node) {
    assert(this->mapping.begin() > this->path_graph.max_node_id());
}

void PhaseUnfolder::unfold(MutableHandleGraph& graph, bool show_progress) {
    
    std::list<bdsg::HashGraph> components = this->complement_components(graph, show_progress);
    
    size_t haplotype_paths = 0;
    bdsg::HashGraph unfolded;
    for (MutableHandleGraph& component : components) {
        haplotype_paths += this->unfold_component(component, graph, unfolded);
    }
    if (show_progress) {
        std::cerr << "Unfolded graph: "
                  << unfolded.get_node_count() << " nodes, " << unfolded.get_edge_count() << " edges on "
                  << haplotype_paths << " paths" << std::endl;
    }
    
    handlealgs::extend(&unfolded, &graph);
}

void PhaseUnfolder::restore_paths(MutableHandleGraph& graph, bool show_progress) const {
    // we include generic to also pick up transcript paths
    this->path_graph.for_each_path_matching({PathSense::GENERIC, PathSense::REFERENCE}, {}, {},
                                            [&](const path_handle_t& path) {
        handle_t prev;
        bool first = true;
        this->path_graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
            handle_t handle = this->path_graph.get_handle_of_step(step);
            vg::id_t id = this->path_graph.get_id(handle);
            handle_t curr;
            if (!graph.has_node(id)) {
                curr = graph.create_handle(this->path_graph.get_sequence(this->path_graph.forward(handle)), id);
                if (this->path_graph.get_is_reverse(handle)) {
                    curr = graph.flip(curr);
                }
            }
            else {
                curr = graph.get_handle(id, this->path_graph.get_is_reverse(handle));
            }
            if (first) {
                // nothing to on the first step
                first = false;
            } else {
                edge_t candidate = make_pair(prev, curr);
                if (!graph.has_edge(candidate)) {
                    graph.create_edge(candidate);
                }
            }
            prev = curr;
        });
    });

    if (show_progress) {
        std::cerr << "Restored graph: " << graph.get_node_count() << " nodes" << std::endl;
    }
}

vg::id_t path_node(const vector<pair<vg::id_t, bool>>& path, size_t i) {
    return path[i].first;
}

vg::id_t path_node(const gbwt::vector_type& path, size_t i) {
    return gbwt::Node::id(path[i]);
}

size_t path_size(const vector<pair<vg::id_t, bool>>& path) {
    return path.size();
}

size_t path_size(const gbwt::vector_type& path) {
    return path.size();
}

bool path_reverse(const vector<pair<vg::id_t, bool>>& path, size_t i) {
    return path[i].second;
}

bool path_reverse(const gbwt::vector_type& path, size_t i) {
    return gbwt::Node::is_reverse(path[i]);
}

struct PathBranch {
    size_t offset;
    size_t curr;    // Branch to choose at offset.
    size_t next;    // Branch to choose at offset + 1.

    void advance() {
        offset++;
        curr = next;
        next = 0;
    }
};

std::ostream& operator<<(std::ostream& out, PathBranch branch) {
    out << "(" << branch.offset << ", " << branch.curr << ", " << branch.next << ")";
    return out;
}

template<class PathType>
bool verify_path(const PathType& path, MutableHandleGraph& unfolded, const hash_map<vg::id_t, std::vector<vg::id_t>>& reverse_mapping) {

    if (path_size(path) < 2) {
        return true;
    }

    // For each branching point: (offset, next duplicate at offset, next duplicate at offset + 1).
    // Initialize with all duplicates of the first node.
    std::vector<PathBranch> branches;
    {
        vg::id_t node_id = path_node(path, 0);
        auto iter = reverse_mapping.find(node_id);
        if (iter != reverse_mapping.end()) {
            for (size_t i = 0; i < iter->second.size(); i++) {
                branches.push_back({ 0, i, 0 });
            }
        } else {
            branches.push_back({ 0, 0, 0 });
        }
    }

    // Note that we can discard all unexplored branches every time the graph
    // contains the original node or only one duplicate of the current node.
    while (!branches.empty()) {
        PathBranch branch = branches.back(); branches.pop_back();
        vg::id_t node_id = path_node(path, branch.offset);
        auto iter = reverse_mapping.find(node_id);
        if (iter != reverse_mapping.end()) {
            node_id = iter->second[branch.curr];
        }
        gbwt::node_type curr = gbwt::Node::encode(node_id, path_reverse(path, branch.offset));

        // Extend the next path from the current branch.
        while (branch.offset + 1 < path_size(path)) {
            node_id = path_node(path, branch.offset + 1);
            size_t duplicates = 0;
            iter = reverse_mapping.find(node_id);
            if (iter != reverse_mapping.end()) {
                node_id = iter->second[branch.next];
                duplicates = iter->second.size();
                if (branch.next + 1 < duplicates) {
                    branches.push_back({ branch.offset, branch.curr, branch.next + 1 });
                }
            }
            gbwt::node_type next = gbwt::Node::encode(node_id, path_reverse(path, branch.offset + 1));
            edge_t candidate = PhaseUnfolder::make_edge(unfolded, curr, next);
            if (!unfolded.has_edge(candidate)) {
                break;
            }
            if (duplicates <= 1) {  // All paths corresponding to this must go through the next node.
                branches.clear();
            }
            curr = next;
            branch.advance();
        }
        if (branch.offset + 1 >= path_size(path)) {
            return true;
        }
    }

    return false;
}

template<class Decoder>
void printId(vg::id_t id) {
    std::cerr << Decoder::id(id);
    if (Decoder::is_reverse(id)) {
        std::cerr << " (reverse)";
    }
}

size_t PhaseUnfolder::verify_paths(MutableHandleGraph& unfolded, bool show_progress) const {

    // Create a mapping from original -> duplicates.
    hash_map<vg::id_t, std::vector<vg::id_t>> reverse_mapping;
    for (gcsa::size_type duplicate = this->mapping.begin(); duplicate < this->mapping.end(); duplicate++) {
        vg::id_t original_id = this->mapping(duplicate);
        reverse_mapping[original_id].push_back(duplicate);
        if (unfolded.has_node(original_id)) {
            reverse_mapping[original_id].push_back(original_id);
        }
    }
    for (auto& mapping : reverse_mapping) {
        gcsa::removeDuplicates(mapping.second, false);
    }

    size_t total_paths = this->path_graph.get_path_count() + this->gbwt_index.sequences(), verified = 0, failures = 0;
    std::set<vg::id_t> failed_threads;
    ProgressBar* progress = nullptr;
    size_t progress_step = std::max(total_paths / 100, static_cast<size_t>(32));
    if (show_progress) {
        progress = new ProgressBar(total_paths, "Verifying paths");
        progress->Progressed(verified);
    }

    this->path_graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        vector<pair<vg::id_t, bool>> path;
        this->path_graph.for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
            handle_t handle = this->path_graph.get_handle_of_step(step);
            path.push_back(make_pair(this->path_graph.get_id(handle),
                                     this->path_graph.get_is_reverse(handle)));
        });
        bool successful = verify_path(path, unfolded, reverse_mapping);
        if (!successful) {
            failures++;
        }
        verified++;
        if (show_progress && (verified % progress_step == 0 || verified >= total_paths)) {
            progress->Progressed(verified);
        }
    });
    
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < this->gbwt_index.sequences(); i++) {

        path_type path = this->gbwt_index.extract(i);
        bool successful = verify_path(path, unfolded, reverse_mapping);
        
#pragma omp critical
        {
            if (!successful) {
                failures++;
                failed_threads.insert(i);
            }
            verified++;
            if (show_progress && (verified % progress_step == 0 || verified >= total_paths)) {
                progress->Progressed(verified);
            }
        }
    }

    if (show_progress) {
        delete progress; progress = nullptr;
        std::cerr << std::endl;

        for (vg::id_t thread_id : failed_threads) {
            path_type path = this->gbwt_index.extract(thread_id);
            std::cerr << "Failed: "; printId<gbwt::Path>(thread_id);
            std::cerr << ": from "; printId<gbwt::Node>(path.front());
            std::cerr << ", to "; printId<gbwt::Node>(path.back());
            std::cerr << ", length " << path.size() << std::endl;
        }
    }

    return failures;
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
    assert(this->mapping.begin() > this->path_graph.max_node_id());
}

vg::id_t PhaseUnfolder::get_mapping(vg::id_t node) const {
    return this->mapping(node);
}

std::list<bdsg::HashGraph> PhaseUnfolder::complement_components(MutableHandleGraph& graph, bool show_progress) {
    
    bdsg::HashGraph complement;

    // checks whether the graph contains an edge
    auto graph_has_edge = [&](const vg::id_t from_id, const vg::id_t to_id,
                              const bool from_rev, const bool to_rev) {
        if (graph.has_node(from_id) && graph.has_node(to_id)) {
            return graph.has_edge(graph.get_handle(from_id, from_rev), graph.get_handle(to_id, to_rev));
        }
        return false;
    };
    
    // checks whether an edge from the PathHandleGraph is in the graph
    auto graph_has_path_graph_edge = [&](const handle_t& from, const handle_t& to) {
        return graph_has_edge(path_graph.get_id(from), path_graph.get_id(to),
                              path_graph.get_is_reverse(from), path_graph.get_is_reverse(to));
    };
    
    // checks whether an edge from the GBWT is in the graph
    auto graph_has_gbwt_edge = [&](const gbwt::node_type& from, const gbwt::node_type& to) {
        return graph_has_edge(gbwt::Node::id(from), gbwt::Node::id(to),
                              gbwt::Node::is_reverse(from), gbwt::Node::is_reverse(to));
    };
    
    // takes a handle to the PathHandleGraph and returns the equivalent handle in
    // the complement, making the node if necessary
    auto get_or_make_complement_handle = [&](const handle_t& counterpart) {
        vg::id_t id = path_graph.get_id(counterpart);
        if (!complement.has_node(id)) {
            complement.create_handle(path_graph.get_sequence(path_graph.forward(counterpart)), id);
        }
        return complement.get_handle(id, path_graph.get_is_reverse(counterpart));
    };
    
    // takes an edge in the XG and ensures that it exists in the complement
    auto make_complement_edge = [&](const handle_t& from, const handle_t& to) {
        handle_t comp_from = get_or_make_complement_handle(from);
        handle_t comp_to = get_or_make_complement_handle(to);
        if (!complement.has_edge(comp_from, comp_to)) {
            complement.create_edge(comp_from, comp_to);
        }
    };
    
    // Add missing edges supported by XG paths.
    this->path_graph.for_each_path_handle([&](const path_handle_t& path) {
        handle_t prev;
        bool first = true;
        this->path_graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
            handle_t handle = this->path_graph.get_handle_of_step(step);
            if (!first) {
                if (!graph_has_path_graph_edge(prev, handle)) {
                    make_complement_edge(prev, handle);
                }
            }
            else {
                first = false;
            }
            prev = handle;
        });
    });

    // Add missing edges supported by GBWT threads, but only if the nodes exist
    // in the original graph.
    for (gbwt::comp_type comp = 1; comp < this->gbwt_index.effective(); comp++) {
        gbwt::node_type gbwt_node = this->gbwt_index.toNode(comp);
        if (!this->path_graph.has_node(gbwt::Node::id(gbwt_node))) {
            continue;
        }
        
        std::vector<gbwt::edge_type> outgoing = this->gbwt_index.edges(gbwt_node);
        for (gbwt::edge_type outedge : outgoing) {
            if (outedge.first == gbwt::ENDMARKER || !this->path_graph.has_node(gbwt::Node::id(outedge.first))) {
                continue;
            }
            if (!graph_has_gbwt_edge(gbwt_node, outedge.first)) {
                make_complement_edge(path_graph.get_handle(gbwt::Node::id(gbwt_node),
                                                           gbwt::Node::is_reverse(gbwt_node)),
                                     path_graph.get_handle(gbwt::Node::id(outedge.first),
                                                           gbwt::Node::is_reverse(outedge.first)));
            }
        }
    }

    std::list<bdsg::HashGraph> components = algorithms::disjoint_components(complement);
    if (show_progress) {
        std::cerr << "Complement graph: "
                  << complement.get_node_count() << " nodes, " << complement.get_edge_count() << " edges in "
                  << components.size() << " components" << std::endl;
    }
    return components;
}

size_t PhaseUnfolder::unfold_component(MutableHandleGraph& component, MutableHandleGraph& graph, MutableHandleGraph& unfolded) {
    // Find the border nodes shared between the component and the graph.
    component.for_each_handle([&](const handle_t& handle) {
        vg::id_t id = component.get_id(handle);
        if (graph.has_node(id)) {
            this->border.insert(id);
        }
    });

    // Generate the paths starting from each border node.
    for (vg::id_t start_node : this->border) {
        this->generate_paths(component, start_node);
    }

    // Generate the threads for each node.
    component.for_each_handle([&](const handle_t& handle) {
        this->generate_threads(component, component.get_id(handle));
    });
    
    auto insert_node = [&](gbwt::node_type node) {
        // create a new node
        if (!unfolded.has_node(gbwt::Node::id(node))) {
            handle_t temp = this->path_graph.get_handle(this->get_mapping(gbwt::Node::id(node)));;
            unfolded.create_handle(this->path_graph.get_sequence(temp), gbwt::Node::id(node));
        }
    };

    // Create the unfolded component from the tries.
    for (auto mapping : this->prefixes) {
        gbwt::node_type from = mapping.first.first, to = mapping.second;
        if (from != gbwt::ENDMARKER) {
            insert_node(from);
        }
        insert_node(to);
        if (from != gbwt::ENDMARKER) {
            unfolded.create_edge(make_edge(unfolded, from, to));
        }
    }
    for (auto mapping : this->suffixes) {
        gbwt::node_type from = mapping.second, to = mapping.first.second;
        insert_node(from);
        if (to != gbwt::ENDMARKER) {
            insert_node(to);
            unfolded.create_edge(make_edge(unfolded, from, to));
        }
    }
    for (auto edge : this->crossing_edges) {
        insert_node(edge.first);
        insert_node(edge.second);
        unfolded.create_edge(make_edge(unfolded, edge.first, edge.second));
    }

    size_t haplotype_paths = this->crossing_edges.size();
    this->border.clear();
    this->reference_paths.clear();
    this->prefixes.clear();
    this->suffixes.clear();
    this->crossing_edges.clear();
    return haplotype_paths;
}

void PhaseUnfolder::generate_paths(MutableHandleGraph& component, vg::id_t from) {

    handle_t from_handle = this->path_graph.get_handle(from);
    this->path_graph.for_each_step_on_handle(from_handle, [&](const step_handle_t& _step) {
        // Forward.
        {
            step_handle_t step = _step;
            handle_t handle = this->path_graph.get_handle_of_step(step);
            vg::id_t id = this->path_graph.get_id(handle);
            bool is_rev = this->path_graph.get_is_reverse(handle);
            gbwt::node_type prev = gbwt::Node::encode(id, is_rev);
            path_type buffer(1, prev);
            while (this->path_graph.has_next_step(step)) {
                step = this->path_graph.get_next_step(step);
                handle = this->path_graph.get_handle_of_step(step);
                id = this->path_graph.get_id(handle);
                is_rev = this->path_graph.get_is_reverse(handle);
                if (!component.has_node(id)) {
                    break;  // Found a maximal path, no matching node.
                }
                gbwt::node_type curr = gbwt::Node::encode(id, is_rev);
                edge_t candidate = make_edge(component, prev, curr);
                if (!component.has_edge(candidate)) {
                    break;  // Found a maximal path, no matching edge.
                }
                buffer.push_back(curr);
                if (this->border.find(gbwt::Node::id(curr)) != this->border.end()) {
                    break;  // Found a border-to-border path.
                }
                prev = curr;
            }
            
            bool to_border = (this->border.find(gbwt::Node::id(buffer.back())) != this->border.end());
            this->reference_paths.push_back(buffer);
            this->insert_path(buffer, true, to_border);
        }

        // Backward.
        {
            step_handle_t step = _step;
            handle_t handle = this->path_graph.get_handle_of_step(step);
            vg::id_t id = this->path_graph.get_id(handle);
            bool is_rev = this->path_graph.get_is_reverse(handle);
            gbwt::node_type prev = gbwt::Node::encode(id, !is_rev);
            path_type buffer(1, prev);
            while (this->path_graph.has_previous_step(step)) {
                step = this->path_graph.get_previous_step(step);
                handle = this->path_graph.get_handle_of_step(step);
                id = this->path_graph.get_id(handle);
                is_rev = this->path_graph.get_is_reverse(handle);
                if (!component.has_node(id)) {
                    break;  // Found a maximal path, no matching node.
                }
                gbwt::node_type curr = gbwt::Node::encode(id, !is_rev);
                edge_t candidate = make_edge(component, prev, curr);
                if (!component.has_edge(candidate)) {
                    break;  // Found a maximal path, no matching edge.
                }
                buffer.push_back(curr);
                if (this->border.find(gbwt::Node::id(curr)) != this->border.end()) {
                    break;  // Found a border-to-border path.
                }
                prev = curr;
            }
            
            bool to_border = (this->border.find(gbwt::Node::id(buffer.back())) != this->border.end());
            this->reference_paths.push_back(buffer);
            this->insert_path(buffer, true, to_border);
        }

    });
}

void PhaseUnfolder::generate_threads(MutableHandleGraph& component, vg::id_t from) {

    bool is_internal = (this->border.find(from) == this->border.end());
    this->create_state(from, false, is_internal);
    this->create_state(from, true, is_internal);

    while (!this->states.empty()) {
        state_type state = this->states.top(); this->states.pop();
        vg::id_t node = gbwt::Node::id(state.first.node);
        bool is_reverse = gbwt::Node::is_reverse(state.first.node);

        if (state.second.size() >= 2 && this->border.find(node) != this->border.end()) {
            if (!is_internal) {
                this->extend_path(state.second);
            }
            continue;   // The path reached a border.
        }

        bool was_extended = false;
        handle_t from = component.get_handle(node, is_reverse);
        component.follow_edges(from, false, [&](const handle_t& handle) {
                was_extended |= this->extend_state(state, component.get_id(handle), component.get_is_reverse(handle));
            });
        component.follow_edges(from, true, [&](const handle_t& handle) {
                was_extended |= this->extend_state(state, component.get_id(handle), !component.get_is_reverse(handle));
            });
        if (!was_extended) {
            this->extend_path(state.second);    // Maximal path.
        }
    }
}

void PhaseUnfolder::create_state(vg::id_t node, bool is_reverse, bool starting) {
    gbwt::node_type gbwt_node = gbwt::Node::encode(node, is_reverse);
    search_type search = (starting ? this->gbwt_index.prefix(gbwt_node) : this->gbwt_index.find(gbwt_node));
    if (search.empty()) {
        return;
    }
    this->states.push(std::make_pair(search, path_type(1, search.node)));
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

PhaseUnfolder::path_type canonical_orientation(const PhaseUnfolder::path_type& path, bool& from_border, bool& to_border) {
    PhaseUnfolder::path_type reverse_complement(path.size(), 0);
    for (size_t i = 0; i < path.size(); i++) {
        reverse_complement[path.size() - 1 - i] = gbwt::Node::reverse(path[i]);
    }
    if (reverse_complement < path) {
        std::swap(from_border, to_border);
        return reverse_complement;
    }
    return path;
}

void PhaseUnfolder::extend_path(const path_type& path) {
    if (path.size() < 2) {
        return;
    }
    bool from_border = (this->border.find(gbwt::Node::id(path.front())) != this->border.end());
    bool to_border = (this->border.find(gbwt::Node::id(path.back())) != this->border.end());
    if (from_border && to_border) {
        this->insert_path(path, from_border, to_border);
        return;
    }

    // We must ensure that the path and its reverse complement are extended in
    // the same way.
    path_type to_extend = canonical_orientation(path, from_border, to_border);

    // Try adding a prefix of a reference path to the beginning of the path.
    // Note that the reverse complement of a reference path is also a
    // reference path.
    if (!from_border) {
        for (size_t ref = 0; ref < this->reference_paths.size(); ref++) {
            const path_type& reference = this->reference_paths[ref];
            bool found = false;
            for (size_t i = 0; i < reference.size(); i++) {
                edge_t candidate = make_edge(path_graph, reference[i], to_extend.front());
                if (this->path_graph.has_edge(candidate.first, candidate.second)) {
                    to_extend.insert(to_extend.begin(), reference.begin(), reference.begin() + i + 1);
                    from_border = true;
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
        }
    }

    // Try adding a suffix of a reference path to the end of the path.
    if (!to_border) {
        for (size_t ref = 0; ref < this->reference_paths.size(); ref++) {
            const path_type& reference = this->reference_paths[ref];
            bool found = false;
            for (size_t i = 0; i < reference.size(); i++) {
                edge_t candidate = make_edge(path_graph, to_extend.back(), reference[i]);
                if (this->path_graph.has_edge(candidate.first, candidate.second)) {
                    to_extend.insert(to_extend.end(), reference.begin() + i, reference.end());
                    to_border = true;
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
        }
    }

    this->insert_path(to_extend, from_border, to_border);
}

void PhaseUnfolder::insert_path(const path_type& path, bool from_border, bool to_border) {

    if (path.size() < 2) {
        return;
    }
    path_type to_insert = canonical_orientation(path, from_border, to_border);

    /*
      Break the path in half. For each prefix / suffix, check if we already
      have a mapping for the next node. If not, create a new duplicate of
      the node and insert the mapping into the corresponding trie. Finally
      create a crossing edge between the full prefix and the full suffix.
    */
    
    // Prefixes.
    gbwt::node_type from = to_insert.front();
    if (!from_border) {
        from = this->get_prefix(gbwt::ENDMARKER, from);
    }
    for (size_t i = 1; i < (to_insert.size() + 1) / 2; i++) {
        from = this->get_prefix(from, to_insert[i]);
    }

    // Suffixes.
    gbwt::node_type to = to_insert.back();
    if (!to_border) {
        to = this->get_suffix(to, gbwt::ENDMARKER);
    }
    for (size_t i = to_insert.size() - 2; i >= (to_insert.size() + 1) / 2; i--) {
        to = this->get_suffix(to_insert[i], to);
    }

    // Crossing edge.
    this->crossing_edges.insert(std::make_pair(from, to));
}


gbwt::node_type PhaseUnfolder::get_prefix(gbwt::node_type from, gbwt::node_type node) {
    std::pair<gbwt::node_type, gbwt::node_type> key(from, node);
    if (this->prefixes.find(key) == this->prefixes.end()) {
        gbwt::size_type new_id = this->mapping.insert(gbwt::Node::id(node));
        this->prefixes[key] = gbwt::Node::encode(new_id, gbwt::Node::is_reverse(node));
    }
    return this->prefixes[key];
}

gbwt::node_type PhaseUnfolder::get_suffix(gbwt::node_type node, gbwt::node_type to) {
    std::pair<gbwt::node_type, gbwt::node_type> key(node, to);
    if (this->suffixes.find(key) == this->suffixes.end()) {
        gbwt::size_type new_id = this->mapping.insert(gbwt::Node::id(node));
        this->suffixes[key] = gbwt::Node::encode(new_id, gbwt::Node::is_reverse(node));
    }
    return this->suffixes[key];
}

} 
