/**
 * \file
 * Implementations for algorithm explanation utilities.
 */

#include "explainer.hpp"

#include <structures/union_find.hpp>

#include <sstream>

namespace vg {

std::atomic<size_t> Explainer::next_explanation_number {0};

bool Explainer::save_explanations = false;

Explainer::Explainer() : explanation_number(Explainer::next_explanation_number++) {
    // Nothing to do!
}

Explainer::~Explainer() {
    // Nothing to do!
}

ProblemDumpExplainer::ProblemDumpExplainer(const std::string& name) : Explainer() {
    if (!Explainer::save_explanations) {
        return;
    }
    out.open(name + std::to_string(explanation_number) + ".json");
}

ProblemDumpExplainer::~ProblemDumpExplainer() {
    // Nothing to do!
}

void ProblemDumpExplainer::object_start() {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "{";
}

void ProblemDumpExplainer::object_end() {
    if (!Explainer::save_explanations) {
        return;
    }
    out << "}";
    need_comma = true;
}

void ProblemDumpExplainer::array_start() {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "[";
}

void ProblemDumpExplainer::array_end() {
    if (!Explainer::save_explanations) {
        return;
    }
    out << "]";
    need_comma = true;
}

void ProblemDumpExplainer::key(const std::string& k) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "\"" << k << "\":";
}

void ProblemDumpExplainer::value(const std::string& v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "\"" << v << "\"";
    need_comma = true;
}

void ProblemDumpExplainer::value(double v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << v;
    need_comma = true;
}

void ProblemDumpExplainer::value(size_t v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "\"" << v << "\"";
    need_comma = true;
}

void ProblemDumpExplainer::value(int v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "\"" << v << "\"";
    need_comma = true;
}

void ProblemDumpExplainer::value(bool v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << (v ? "true" : "false");
    need_comma = true;
}

void ProblemDumpExplainer::value(vg::id_t v) {
    if (!Explainer::save_explanations) {
        return;
    }
    comma();
    out << "\"" << v << "\"";
    need_comma = true;
}

void ProblemDumpExplainer::value(const pos_t& v) {
    if (!Explainer::save_explanations) {
        return;
    }
    object_start();
    key("node_id");
    value(id(v));
    if (offset(v) != 0) {
        key("offset");
        value(offset(v));
    }
    if (is_rev(v)) {
        key("is_reverse");
        value(true);
    }
    object_end();
}

void ProblemDumpExplainer::value(const HandleGraph& v) {
    if (!Explainer::save_explanations) {
        return;
    }
    object_start();
    key("node");
    array_start();
    v.for_each_handle([&](const handle_t& h) {
        // Put all the nodes in the node array
        object_start();
        key("id");
        value(v.get_id(h));
        key("sequence");
        value(v.get_sequence(h));
        object_end();
    });
    array_end();
    key("edge");
    array_start();
    v.for_each_edge([&](const edge_t& e) {
        // Put all the edges in the edge array
        object_start();
        key("from");
        value(v.get_id(e.first));
        if (v.get_is_reverse(e.first)) {
            key("from_start");
            value(true);
        }
        key("to");
        value(v.get_id(e.second));
        if (v.get_is_reverse(e.second)) {
            key("to_end");
            value(true);
        }
        object_end();
    });
    array_end();
    object_end();
}

void ProblemDumpExplainer::value(const handle_t& v, const HandleGraph& context) {
    if (!Explainer::save_explanations) {
        return;
    }
    // Implement via pos_t serialization.
    this->value(make_pos_t(context.get_id(v), context.get_is_reverse(v), 0));
}

const size_t DiagramExplainer::MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY {5};

DiagramExplainer::DiagramExplainer() : Explainer() {
    // Nothing to do!
}

DiagramExplainer::~DiagramExplainer() {
    if (!Explainer::save_explanations) {
        return;
    }
    write_connected_components();
}

void DiagramExplainer::add_globals(const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }
    std::copy(annotations.begin(), annotations.end(), std::back_inserter(globals));
}

void DiagramExplainer::add_node(const std::string& id, const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }
    nodes.emplace(id, annotations);
}

void DiagramExplainer::ensure_node(const std::string& id, const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }
    auto found = nodes.find(id);
    if (found == nodes.end()) {
        nodes.emplace_hint(found, id, annotations);
    }
}

void DiagramExplainer::add_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }
    edges.emplace(std::make_pair(a_id, b_id), annotations);
}

void DiagramExplainer::ensure_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }
    auto key = std::make_pair(a_id, b_id);
    auto found = edges.find(key);
    if (found == edges.end()) {
        edges.emplace_hint(found, std::move(key), annotations);
    }
}

void DiagramExplainer::suggest_edge(const std::string& a_id, const std::string& b_id, const std::string& category, double importance, const annotation_t& annotations) {
    if (!Explainer::save_explanations) {
        return;
    }

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

void DiagramExplainer::for_each_edge(const std::function<void(const edge_ref_t&)>& iteratee) const {
    for (auto& kv : edges) {
        // Do all the required edges
        iteratee(edge_ref_t(kv.first.first, kv.first.second, kv.second));
    }
    
    for (auto& kv : suggested_edges) {
        for (auto& suggestion : kv.second) {
            const stored_edge_t& edge = suggestion.second;
            // Do all the surviving suggested edges
            iteratee(edge_ref_t(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge)));
        }
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
    
    for_each_edge([&](const edge_ref_t& edge) {
        // Connect connected components for each edge
        components.union_groups(id_to_index.at(std::get<0>(edge)), id_to_index.at(std::get<1>(edge)));
    });
    
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
    
    for_each_edge([&](const edge_ref_t& edge) {
        // Add each edge to the file for its group
        size_t group = components.find_group(id_to_index.at(std::get<0>(edge)));
        write_edge(files_by_group.at(group), std::get<0>(edge), std::get<1>(edge), std::get<2>(edge));
    });
    
    for (auto& kv : files_by_group) {
        // Close out all the files
        kv.second << "}" << std::endl;
        kv.second.close();
    }
}

}
