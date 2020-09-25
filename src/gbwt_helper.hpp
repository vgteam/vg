#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWT.
 */

#include <vector>

#include "position.hpp"

#include <gbwt/dynamic_gbwt.h>
#include <handlegraph/mutable_path_handle_graph.hpp>

namespace vg {

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples);

//------------------------------------------------------------------------------

/// Convert gbwt::node_type to handle_t.
inline handle_t gbwt_to_handle(const HandleGraph& graph, gbwt::node_type node) {
    return graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
}

/// Convert gbwt::node_type and an offset as size_t to pos_t.
inline pos_t gbwt_to_pos(gbwt::node_type node, size_t offset) {
    return make_pos_t(gbwt::Node::id(node), gbwt::Node::is_reverse(node), offset);
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

/// Convert Path to a GBWT path.
inline gbwt::vector_type path_to_gbwt(const Path& path) {
    gbwt::vector_type result(path.mapping_size());
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = mapping_to_gbwt(path.mapping(i));
    }
    return result;
}

/// Extract a path as a GBWT path.
gbwt::vector_type extract_as_gbwt_path(const PathHandleGraph& graph, const std::string& path_name);

// Find all predecessor nodes of the path, ignoring self-loops.
gbwt::vector_type path_predecessors(const PathHandleGraph& graph, const std::string& path_name);

//------------------------------------------------------------------------------

// GBWT construction helpers.

/// Determine the node width in bits for the GBWT nodes based on the given graph.
gbwt::size_type gbwt_node_width(const HandleGraph& graph);

/// Finish GBWT construction and optionally print the metadata.
void finish_gbwt_constuction(gbwt::GBWTBuilder& builder,
    const std::vector<std::string>& sample_names,
    const std::vector<std::string>& contig_names,
    size_t haplotype_count, bool print_metadata);

//------------------------------------------------------------------------------

/// Insert a GBWT thread into the graph and return its name. Returns an empty string on failure.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string insert_gbwt_path(MutablePathHandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id);

/// Extract a GBWT thread as a path in the given graph.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
Path extract_gbwt_path(const HandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id);

/// Get a string representation of a thread name stored in GBWT metadata.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string thread_name(const gbwt::GBWT& gbwt_index, gbwt::size_type id);

//------------------------------------------------------------------------------

/// Transform the paths into a GBWT index. Primarily for testing.
gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
