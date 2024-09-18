#ifndef VG_EXPLAINER_HPP_INCLUDED
#define VG_EXPLAINER_HPP_INCLUDED

/**
 * \file
 * Contains utility classes for producing algorithms which can explain
 * themselves and capture monitoring statistics in an efficient way.
 */

#include <atomic>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <functional>
#include <queue>

// For pair hash overload
#include "hash_map.hpp"

#include "types.hpp"
#include "handle.hpp"

namespace vg {

/**
 * Base explainer class. Handles making sure each explanation has a different unique number.
 */
class Explainer {
public:
    /// Determine if explanations should be generated.
    static bool save_explanations;

    /// Construct an Explainer that will save to one or more files
    Explainer(bool enabled);
    
    /// Close out the files being explained to
    virtual ~Explainer();

    /// Conversion to bool so you can use an explainer as a condition on code
    /// to write to it.
    inline operator bool() const {
        return explaining();
    }

protected:
    /// What number explanation are we? Distinguishes different objects.
    size_t explanation_number;
    
    /// Determines if this explainer should generate explanations.
    bool enabled;

    /// Counter used to give different explanations their own unique filenames.
    static std::atomic<size_t> next_explanation_number;

    /// Function to check if we should be explaining.
    inline bool explaining() const {
        return this->enabled && Explainer::save_explanations;
    }
};

/**
 * Widget to log a TSV of data as an explanation.
 */
class TSVExplainer : public Explainer {
public:
    /// Construct a TSVExplainer that will save a table to a file.
    TSVExplainer(bool enabled, const std::string& name = "data");
    /// Close out the file being explained to
    ~TSVExplainer();

    /// Start a new line. Must call this before field().
    void line();

    /// Add a field with a string value
    void field(const std::string& value);

    /// Add a field with an integral value
    void field(size_t value);

protected:
    /// Stream being written to
    ofstream out;
    /// Whether we need a tab befroe the next value
    bool need_tab = false;
    /// Whether we need a newline before the next line
    bool need_line = false;
};

/**
 * Widget to serialize somewhat structured logs.
 */
class ProblemDumpExplainer : public Explainer {
public:
    /// Construct a ProblemDumpExplainer that will save a dump of a problem to a file.
    ProblemDumpExplainer(bool enabled, const std::string& name = "problem");
    /// Close out the file being explained to
    ~ProblemDumpExplainer();
    
    // We think in JSON, but with support for vg types.
    
    /// Begin an object in a value context.
    void object_start();
    /// End an object after its last value.
    void object_end();
    /// Begin an array in a value context.
    void array_start();
    /// End an array after its last value.
    void array_end();
    
    /// Put the key for a value, inside an object
    void key(const std::string& k);
    
    /// Put a value after a key or in an array.
    /// Assumes the string is pre-escaped.
    void value(const std::string& v);
    /// Put a value after a key or in an array.
    void value(double v);
    /// Put a value after a key or in an array.
    void value(size_t v);
    /// Put a value after a key or in an array.
    void value(int v);
    /// Put a value after a key or in an array.
    void value(bool v);
    /// Put a value after a key or in an array.
    void value(vg::id_t v);
    /// Put a value after a key or in an array.
    /// Represents the position as a vg Protobuf Position.
    void value(const pos_t& v);
    /// Put a value after a key or in an array.
    /// Represents the graph as a single chunk vg Protobuf Graph.
    void value(const HandleGraph& v);
    /// Put a value after a key or in an array.
    /// Represents a handle as a vg Protobuf Position.
    void value(const handle_t& v, const HandleGraph& context);
    
protected:
    /// Stream being written to
    ofstream out;
    /// Whether we need a comma before the next key or value.
    bool need_comma = false;
    
    /// Write a separating comma if needed.
    inline void comma() {
        if (need_comma) {
            out << ",";
            need_comma = false;
        }
    }
};

/**
 * Widget to log statistics to a GraphViz file.
 */
class DiagramExplainer : public Explainer {
public:
    // We define a type for miscelaneous annotations, since we don't have kwargs
    using annotation_t = std::vector<std::pair<std::string, std::string>>;

    /// Construct a DiagramExplainer that will save a diagram to one or more files.
    DiagramExplainer(bool enabled);
    /// Close out the files being explained to
    ~DiagramExplainer();
    
    /// Add global annotations (like rankdir)
    void add_globals(const annotation_t& annotations);
    
    /// Add a node. Optionally give it the given annotation, which must be pre-escaped.
    /// The node is assumed not to exist already
    void add_node(const std::string& id, const annotation_t& annotations = {});
    
    /// Add a node. Optionally give it the given annotation, which must be pre-escaped.
    /// Deduplicates multiple calls with the same ID.
    void ensure_node(const std::string& id, const annotation_t& annotations = {});
    
    /// Add an edge. Optionally give it the given annotation, which must be
    /// pre-escaped. The edge is assumed not to exist already.
    void add_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations = {});
    
    /// Add an edge. Optionally give it the given annotation, which must be
    /// pre-escaped. Deduplicates multiple calls with the same IDs in the same order.
    void ensure_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations = {});
    
    /// Add an optional edge. Optionally give it the given annotation, which must be pre-escaped.
    /// Only the k most important edges in each category will actually render
    void suggest_edge(const std::string& a_id, const std::string& b_id, const std::string& category, double importance, const annotation_t& annotations = {});
    
protected:
    /// Collection of all global diagram key-value pairs (like rankdir)
    annotation_t globals;

    /// Collection of all nodes, by ID
    std::unordered_map<std::string, annotation_t> nodes;
    
    /// We will need to store edges
    using stored_edge_t = std::tuple<std::string, std::string, annotation_t>;
    /// And show them to people
    using edge_ref_t = std::tuple<const std::string&, const std::string&, const annotation_t&>;
    
    /// Collection of all required edges
    std::unordered_map<std::pair<std::string, std::string>, annotation_t> edges;
    
    using suggested_edge_t = std::pair<double, stored_edge_t>;
    
    /// Top k most important edges for each suggested edge category
    std::unordered_map<std::string, std::vector<suggested_edge_t>> suggested_edges;
    
    /// Limit on suggested edges
    static const size_t MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY;
    
    /// Loop over all the edges, across all kinds of storage
    void for_each_edge(const std::function<void(const edge_ref_t&)>& iteratee) const;
    
    /// Save the annotations for a node or edge, if any.
    void write_annotations(std::ostream& out, const annotation_t& annotations) const;
    
    /// Write out a node
    void write_node(std::ostream& out, const std::string& id, const annotation_t& annotations) const;
    
    /// Write out an edge
    void write_edge(std::ostream& out, const std::string& a_id, const std::string& b_id, const annotation_t& annotations) const;
    
    /// Write out globals
    void write_globals(std::ostream& out, const annotation_t& annotations) const;
    
    /// Write each connected component to a different file
    void write_connected_components() const;
    
};

/**
 * Explainer that can dump anything that has a:
 *   void to_dot(ostream& out) const;
 * method, such as a Funnel.
 */
template<typename T>
class DotDumpExplainer : public Explainer {
public:
    /// Construct a DotDumpExplainer that will save a diagram to a file
    DotDumpExplainer(bool enabled, const T& to_dump);
};

template<typename T>
DotDumpExplainer<T>::DotDumpExplainer(bool enabled, const T& to_dump) : Explainer(enabled) {
    if (!explaining()) {
        return;
    }
    // Open the dot file
    std::ofstream out("dotdump" + std::to_string(explanation_number) + ".dot");
    // And dump to it
    to_dump.to_dot(out);
}

/**
 * Explainer that can dump a handle graph.
 */
class SubgraphExplainer: public Explainer {
public:

    /// Construct an explainer that will save a single graph.
    SubgraphExplainer(bool enabled);

    /// Write out a subgraph.
    void subgraph(const HandleGraph& graph);
};


}
 
#endif
