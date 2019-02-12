#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility functions for working with GBWT.
 */

#include "handle.hpp"
#include "position.hpp"
#include "xg.hpp"

#include <gbwt/dynamic_gbwt.h>

namespace vg {

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
    return { (int64_t)(gbwt::Node::id(node)), gbwt::Node::is_reverse(node) };
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

/// Traverse all haplotype-consistent kmers in the graph and call lambda() for each kmer.
/// Uses multiple threads, so the lambda should be thread-safe.
void for_each_kmer(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t k,
                   const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                   bool parallel);

/// Traverse all haplotype-consistent window in the graph and call lambda() for each kmer.
/// Uses multiple threads, so the lambda should be thread-safe.
/// A window starts with the sequence of a node and is followed by window_size - 1 bases
/// from subsequent nodes. If no extensions are possible, a shorter substring of
/// length >= window_size also qualifies as a window.
void for_each_window(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t window_size,
                     const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                     bool parallel);

/// Iterate over all windows in the graph, running lambda on each.
void for_each_window(const HandleGraph& graph, size_t window_size,
                     const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                     bool parallel);

/// Transform the paths into a GBWT index. Primarily for testing.
gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths);

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
