#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility functions for working with GBWT.
 */

#include "handle.hpp"
#include "position.hpp"

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

/// Convert handle_t to gbwt::node_type.
inline gbwt::node_type handle_to_gbwt(const HandleGraph& graph, handle_t handle) {
    return gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
}

/// Extract gbwt::node_type from pos_t.
inline gbwt::node_type pos_to_gbwt(pos_t pos) {
    return gbwt::Node::encode(std::get<0>(pos), std::get<1>(pos));
}

/// Stores haplotype-consistent traversal of the graph and the corresponding sequence.
struct GBWTTraversal {
    /// The traversal as a sequence of (begin, end) pairs.
    std::vector<std::pair<pos_t, pos_t>> traversal;
    /// The sequence.
    std::string seq;
    /// GBWT search state at the end of the traversal.
    gbwt::SearchState state;
};

/// Traverse all haplotype-consistent kmers in the graph and call lambda() for each kmer.
void for_each_kmer(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t k,
                   const function<void(const GBWTTraversal&)>& lambda);

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
