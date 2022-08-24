/**
 * \file
 * Implementations for algorithm explanation utilities.
 */

#include "explainer.hpp"

namespace vg {

std::atomic<size_t> DiagramExplainer::next_explanation_number {0};
const size_t DiagramExplainer::MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY {5};

DiagramExplainer::DiagramExplainer(const annotation_t& global_annotations) {
    // Get a unique filename
    file = std::ofstream("graph" + std::to_string(DiagramExplainer::next_explanation_number++) + ".dot");
    // And start the graphviz graph
    file << "digraph explanation {" << std::endl;
    for (auto& kv : global_annotations) {
        // Add all the globals like rankdir
        file << kv.first << "=\"" << kv.second << "\";" << std::endl;
    }
}

DiagramExplainer::~DiagramExplainer() {
    // Finish with the best suggested edges for each category
    for (auto& kv : suggested_edges) {
        auto& queue = kv.second;
        while (!queue.empty()) {
            // Unpack all the surviving suggestions for the category
            auto& suggestion = queue.top().second;
            auto& a_id = std::get<0>(suggestion);
            auto& b_id = std::get<1>(suggestion);
            auto& annotations = std::get<2>(suggestion);
            // And turn those into real edges
            add_edge(a_id, b_id, annotations);
            queue.pop();
        }
    }
    // End the graphviz grap
    file << "}" << std::endl;
}

void DiagramExplainer::add_node(const std::string& id, const annotation_t& annotations) {
    file << id;
    write_annotations(annotations);
    file << ";" << std::endl;
}

void DiagramExplainer::ensure_node(const std::string& id, const annotation_t& annotations) {
    auto found = ensured_node_ids.find(id);
    if (found == ensured_node_ids.end()) {
        ensured_node_ids.emplace_hint(found, id);
        add_node(id, annotations);
    }
}

void DiagramExplainer::add_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations) {
    file << a_id << " -> " << b_id;
    write_annotations(annotations);
    file << ";" << std::endl;
}

void DiagramExplainer::suggest_edge(const std::string& a_id, const std::string& b_id, const std::string& category, double importance, const annotation_t& annotations) {
    // Find the queue it goes in
    auto& queue = suggested_edges[category];
    
    // Put the suggestion in
    queue.emplace(importance, std::make_tuple(a_id, b_id, annotations));
    while (queue.size() > DiagramExplainer::MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY) {
        // Clamp to size limit
        queue.pop();
    }
    
}

void DiagramExplainer::write_annotations(const annotation_t& annotations) {
    if (!annotations.empty()) {
        file << " [";
        for (auto& annotation : annotations) {
            file << annotation.first << "=\"" << annotation.second << "\", ";
        }
        file << "]";
    }
}

}
