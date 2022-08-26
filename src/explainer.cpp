/**
 * \file
 * Implementations for algorithm explanation utilities.
 */

#include "explainer.hpp"

#include <structures/union_find.hpp>

#include <sstream>

namespace vg {

std::atomic<size_t> DiagramExplainer::next_explanation_number {0};
const size_t DiagramExplainer::MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY {5};

DiagramExplainer::DiagramExplainer() : explanation_number(DiagramExplainer::next_explanation_number++) {
    // Nothing to do!
}

void DiagramExplainer::add_globals(const annotation_t& annotations) {
    std::copy(annotations.begin(), annotations.end(), std::back_inserter(globals));
}

DiagramExplainer::~DiagramExplainer() {
    write_connected_components();
}

void DiagramExplainer::add_node(const std::string& id, const annotation_t& annotations) {
    nodes.emplace(id, annotations);
}

void DiagramExplainer::ensure_node(const std::string& id, const annotation_t& annotations) {
    auto found = nodes.find(id);
    if (found == nodes.end()) {
        nodes.emplace_hint(found, id, annotations);
    }
}

void DiagramExplainer::add_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations) {
    edges.emplace_back(a_id, b_id, annotations);
}

void DiagramExplainer::suggest_edge(const std::string& a_id, const std::string& b_id, const std::string& category, double importance, const annotation_t& annotations) {
    // Find the heap it goes in
    auto& heap = suggested_edges[category];
    
    // And make a comparator
    std::greater<suggested_edge_t> comp;
    
    // Put the suggestion in
    heap.emplace_back(importance, std::make_tuple(a_id, b_id, annotations));
    std::push_heap(heap.begin(), heap.end(), comp);
    while (heap.size() > DiagramExplainer::MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY) {
        // Clamp to size limit
        std::pop_heap(heap.begin(), heap.end(), comp);
        heap.pop_back();
    }
    
}

void DiagramExplainer::write_annotations(std::ostream& out, const annotation_t& annotations) const {
    if (!annotations.empty()) {
        out << " [";
        for (auto& annotation : annotations) {
            // Add all the annotations to the thing, in the brackets that we need
            out << annotation.first << "=\"" << annotation.second << "\", ";
        }
        out << "]";
    }
}

void DiagramExplainer::write_node(std::ostream& out, const std::string& id, const annotation_t& annotations) const {
    out << id;
    write_annotations(out, annotations);
    out << ";" << std::endl;
}

void DiagramExplainer::write_edge(std::ostream& out, const std::string& a_id, const std::string& b_id, const annotation_t& annotations) const {
    out << a_id << " -> " << b_id;
    write_annotations(out, annotations);
    out << ";" << std::endl;
}

void DiagramExplainer::write_globals(std::ostream& out, const annotation_t& annotations) const {
    for (auto& kv : annotations) {
        // Add all the globals, which have a different syntax than item annotations.
        out << kv.first << "=\"" << kv.second << "\";" << std::endl;
    }
}

void DiagramExplainer::write_connected_components() const {
    // Number all the nodes
    std::vector<decltype(nodes)::const_iterator> node_order;
    // TODO: Use a symbol-registering widget so we don't need to keep both these tables.
    std::unordered_map<std::string, size_t> id_to_index;
    node_order.reserve(nodes.size());
    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        id_to_index.emplace(it->first, node_order.size());
        node_order.push_back(it);
    }
    
    // Compose connected components
    structures::UnionFind components(node_order.size());
    for (const stored_edge_t& edge : edges) {
        // Connect connected components for each required edge
        components.union_groups(id_to_index.at(std::get<0>(edge)), id_to_index.at(std::get<1>(edge)));
    }
    
    for (auto& kv : suggested_edges) {
        for (auto& suggestion : kv.second) {
            const stored_edge_t& edge = suggestion.second;
            // Connect connected components for each surviving suggested edge
            components.union_groups(id_to_index.at(std::get<0>(edge)), id_to_index.at(std::get<1>(edge)));
        }
    }
    
    std::unordered_map<size_t, std::ofstream> files_by_group;
    for (size_t i = 0; i < node_order.size(); i++) {
        // For each node
        
        // Make sure we have a file for the connected component it goes in
        size_t group = components.find_group(i);
        auto file_it = files_by_group.find(group);
        if (file_it == files_by_group.end()) {
            // We need to open and set up a new file
            std::stringstream name_stream;
            name_stream << "graph" << explanation_number << "-" << files_by_group.size() << ".dot";
            file_it = files_by_group.emplace_hint(file_it, group, name_stream.str());
            
            // Start off with the heading
            file_it->second << "digraph explanation {" << std::endl;
            // And any globals
            write_globals(file_it->second, globals);
        }
        
        // Write the node
        write_node(file_it->second, node_order[i]->first, node_order[i]->second);
    }
    
    for (const stored_edge_t& edge : edges) {
        // Add each required edge to the file for its group
        size_t group = components.find_group(id_to_index.at(std::get<0>(edge)));
        write_edge(files_by_group.at(group), std::get<0>(edge), std::get<1>(edge), std::get<2>(edge));
    }
    
    for (auto& kv : suggested_edges) {
        for (auto& suggestion : kv.second) {
            const stored_edge_t& edge = suggestion.second;
            // Add each surviving suggested edge to the file for its group
            size_t group = components.find_group(id_to_index.at(std::get<0>(edge)));
            write_edge(files_by_group.at(group), std::get<0>(edge), std::get<1>(edge), std::get<2>(edge));
        }
    }
    
    for (auto& kv : files_by_group) {
        // Close out all the files
        kv.second << "}" << std::endl;
        kv.second.close();
    }
}

}
