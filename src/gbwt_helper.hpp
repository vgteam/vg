#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWT.
 */

#include <vector>

#include "handle.hpp"
#include "position.hpp"
#include "xg.hpp"

#include <gbwt/dynamic_gbwt.h>
#include <sdsl/int_vector.hpp>

namespace vg {

//------------------------------------------------------------------------------

/// Convert gbwt::node_type to handle_t.
inline handle_t gbwt_to_handle(const HandleGraph& graph, gbwt::node_type node) {
    return graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
}

/// Convert gbwt::node_type and an offset as size_t to pos_t.
inline pos_t gbwt_to_pos(gbwt::node_type node, size_t offset) {
    return make_pos_t(gbwt::Node::id(node), gbwt::Node::is_reverse(node), offset);
}

/// Convert gbwt::node_type to xg::XG:ThreadMapping.
inline xg::XG::ThreadMapping gbwt_to_thread_mapping(gbwt::node_type node) {
    return { static_cast<int64_t>(gbwt::Node::id(node)), gbwt::Node::is_reverse(node) };
}

/// Convert handle_t to gbwt::node_type.
inline gbwt::node_type handle_to_gbwt(const HandleGraph& graph, handle_t handle) {
    return gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
}

/// Extract gbwt::node_type from pos_t.
inline gbwt::node_type pos_to_gbwt(pos_t pos) {
    return gbwt::Node::encode(id(pos), is_rev(pos));
}

/// Convert Mapping to gbwt::node_type.
inline gbwt::node_type mapping_to_gbwt(const Mapping& mapping) {
    return gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse());
}

/// Convert a node on xg::XGPath to gbwt::node_type.
inline gbwt::node_type xg_path_to_gbwt(const xg::XGPath& path, size_t i) {
    return gbwt::Node::encode(path.node(i), path.is_reverse(i));
}

//------------------------------------------------------------------------------

/**
 * A runtime-only HandleGraph implementation that uses GBWT for graph topology and
 * extracts sequences from another HandleGraph. Faster sequence access but slower
 * graph navigation than in XG. Also supports a version of follow_edges() that
 * takes only paths supported by the indexed haplotypes.
 */
class GBWTGraph : public HandleGraph {
public:
    /// Create a graph backed by the GBWT index and extract the sequences from the
    /// given HandleGraph.
    GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source);

    /// Copy constructor.
    GBWTGraph(const GBWTGraph& source);

    /// Move constructor.
    GBWTGraph(GBWTGraph&& source);

    const gbwt::GBWT&   index;
    std::vector<char>   sequences;
    sdsl::int_vector<0> offsets;
    std::vector<bool>   real_nodes;
    size_t              total_nodes;

//------------------------------------------------------------------------------

    // Standard HandleGraph interface.

    // Method to check if a node exists by ID
    virtual bool has_node(id_t node_id) const;

    /// Look up the handle for the node with the given ID in the given orientation.
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;

    // Copy over the visit version which would otherwise be shadowed.
    using HandleGraph::get_handle;

    /// Get the ID from a handle.
    virtual id_t get_id(const handle_t& handle) const;

    /// Get the orientation of a handle.
    virtual bool get_is_reverse(const handle_t& handle) const;

    /// Invert the orientation of a handle (potentially without getting its ID).
    virtual handle_t flip(const handle_t& handle) const;

    /// Get the length of a node.
    virtual size_t get_length(const handle_t& handle) const;

    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual std::string get_sequence(const handle_t& handle) const;

    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns false if we stopped early.
    virtual bool follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;

    // Copy over the template version.
    using HandleGraph::follow_edges;

    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual void for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;

    // Copy over the template version.
    using HandleGraph::for_each_handle;

    /// Return the number of nodes in the graph.
    virtual size_t node_size() const;

    /// Get the minimum node ID used in the graph, if any are used.
    virtual id_t min_node_id() const;

    /// Get the maximum node ID used in the graph, if any are used.
    virtual id_t max_node_id() const;

//------------------------------------------------------------------------------

    // GBWTGraph specific interface.

    /// Convert gbwt::node_type to handle_t.
    static handle_t node_to_handle(gbwt::node_type node) { return as_handle(node); }

    /// Convert handle_t to gbwt::node_type.
    static gbwt::node_type handle_to_node(const handle_t& handle) { return as_integer(handle); }

    /// Convert handle_t to gbwt::SearchState.
    gbwt::SearchState get_state(const handle_t& handle) const { return this->index.find(as_integer(handle)); }

    /// Convert handle_t to gbwt::BidirectionalState.
    gbwt::BidirectionalState get_bd_state(const handle_t& handle) const { return this->index.bdFind(as_integer(handle)); }

    /// Visit all successor states of this state and call iteratee for the state.
    /// Stop and return false if the iteratee returns false.
    /// Note that the state may be empty if no path continues to that node.
    bool follow_edges(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const;

    /// Visit all predecessor/successor states of this state and call iteratee for the state.
    /// Stop and return false if the iteratee returns false.
    /// Note that the state may be empty if no path continues to that node.
    /// Each state corresponds to a path. Going backward extends the path left, while going
    /// extends it right.
    bool follow_edges(gbwt::BidirectionalState state, bool backward,
                      const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const;

private:
    size_t node_offset(gbwt::node_type node) const { return node - this->index.firstNode(); }
    size_t node_offset(const handle_t& handle) const { return this->node_offset(as_integer(handle)); }
    bool is_real(gbwt::node_type node) const { return this->real_nodes[this->node_offset(node) / 2]; }
};

//------------------------------------------------------------------------------

/// Traverse all haplotype-consistent kmers in the graph and call lambda() for each kmer.
/// Uses multiple threads, so the lambda should be thread-safe.
void for_each_kmer(const GBWTGraph& graph, size_t k,
                   const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                   bool parallel);

/// Traverse all haplotype-consistent window in the graph and call lambda() for each kmer.
/// Uses multiple threads, so the lambda should be thread-safe.
/// A window starts with the sequence of a node and is followed by window_size - 1 bases
/// from subsequent nodes. If no extensions are possible, a shorter substring of
/// length >= window_size also qualifies as a window.
void for_each_haplotype_window(const GBWTGraph& graph, size_t window_size,
                               const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                               bool parallel);

/// Iterate over all windows in the graph, running lambda on each.
void for_each_window(const HandleGraph& graph, size_t window_size,
                     const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                     bool parallel);

/// Transform the paths into a GBWT index. Primarily for testing.
gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
