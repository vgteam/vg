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

namespace vg {

/**
 * Widget to log statistics to a GraphViz file.
 */
class DiagramExplainer {
public:
    // We define a type for miscelaneous annotations, since we don't have kwargs
    using annotation_t = std::vector<std::pair<std::string, std::string>>;

    /// Construct a DiagramExplainer that will save a diagram to a file
    DiagramExplainer();
    /// Close out the file being explained to
    ~DiagramExplainer();
    
    /// Add global annotations (like rankdir)
    void add_globals(const annotation_t& annotations);
    
    /// Add a node. Optionally give it the given annotation, which must be pre-escaped.
    /// The node is assumed not to exist already
    void add_node(const std::string& id, const annotation_t& annotations = {});
    
    /// Add a node. Optionally give it the given annotation, which must be pre-escaped.
    /// Deduplicates multiple calls with the same ID.
    void ensure_node(const std::string& id, const annotation_t& annotations = {});
    
    /// Add an edge. Optionally give it the given annotation, which must be pre-escaped.
    void add_edge(const std::string& a_id, const std::string& b_id, const annotation_t& annotations = {});
    
    /// Add an optional edge. Optionally give it the given annotation, which must be pre-escaped.
    /// Only the k most important edges in each category will actually render
    void suggest_edge(const std::string& a_id, const std::string& b_id, const std::string& category, double importance, const annotation_t& annotations = {});
    
protected:
    /// The file we are reporting to
    std::ofstream file;
    
    /// Set used for deduplicating nodes for ensure_node
    std::unordered_set<std::string> ensured_node_ids;
    
    // We will need to store edges for managing suggested edges
    using stored_edge_t = std::tuple<std::string, std::string, annotation_t>;
    using suggested_edge_t = std::pair<double, stored_edge_t>;
    
    /// Top k most important edges for each suggested edge category
    std::unordered_map<std::string, std::priority_queue<suggested_edge_t, std::vector<suggested_edge_t>, std::greater<suggested_edge_t>>> suggested_edges;
    
    
    /// Counter used to give different explanations their own unique filenames.
    static std::atomic<size_t> next_explanation_number;
    
    /// Limit on suggested edges
    static const size_t MAX_DISPLAYED_SUGGESTIONS_PER_CATEGORY;
    
    
    /// Save the annotations for a node or edge, if any.
    void write_annotations(const annotation_t& annotations);
};

}
 
#endif
