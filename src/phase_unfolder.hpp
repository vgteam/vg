#ifndef VG_PHASE_UNFOLDER_HPP_INCLUDED
#define VG_PHASE_UNFOLDER_HPP_INCLUDED

/** \file
 * This file contains the PhaseUnfolder, which replaces pruned regions of
 * the graph with paths supported by XG paths or GBWT threads.
 */

#include "vg.hpp"
#include "handle.hpp"
#include "hash_map.hpp"
#include "gbwt_helper.hpp"

#include <algorithm>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include <gcsa/support.h>
#include <bdsg/hash_graph.hpp>

namespace vg {

/**
 * Transforms the pruned subregions of the input graph into collections of
 * disconnected distinct traversal haplotypes. Use in combination with
 * pruning to simplify the graph for GCSA2 indexing without losing observed
 * variation.
 * Requires the XG index of the original graph and an empty GBWT index or
 * an GBWT index of the original graph.
 * Note: PhaseUnfolder only considers paths of length >= 2.
 */
class PhaseUnfolder {
public:
    typedef gbwt::SearchState                 search_type;
    typedef gbwt::vector_type                 path_type;
    typedef std::pair<search_type, path_type> state_type;

    /**
     * Make a new PhaseUnfolder backed by the given XG and GBWT indexes.
     * These indexes must represent the same original graph. 'next_node' should
     * usually be max_node_id() + 1 in the original graph.
     */
    PhaseUnfolder(const PathHandleGraph& path_graph, const gbwt::GBWT& gbwt_index, vg::id_t next_node);

    /**
     * Unfold the pruned regions in the input graph:
     *
     * - Determine the connected components of edges missing from the input
     * graph, as determined by the XG paths and GBWT threads.
     *
     * - For each component, find all border-to-border paths and threads
     * supported by the indexes. Then unfold the component by duplicating the
     * nodes, so that the paths are disjoint, except for their shared prefixes
     * and suffixes.
     *
     * - Extend the input graph with the unfolded components.
     */
    void unfold(MutableHandleGraph& graph, bool show_progress = false);

    /**
     * Restore the edges on XG paths. This is effectively the same as
     * unfolding with an empty GBWT index, except that the inserted nodes will
     * have their original identifiers.
     */
    void restore_paths(MutableHandleGraph& graph, bool show_progress = false) const;

    /**
     * Verify that the graph contains the XG paths and the GBWT threads in the
     * backing indexes. Returns the number of paths for which the verification
     * failed. Uses OMP threads.
     */
    size_t verify_paths(MutableHandleGraph& unfolded, bool show_progress = false) const;

    /**
     * Write the mapping to the specified file with a header. The file will
     * contain mappings from header.next_node - header.mapping_size
     * (inclusive) to header.next_node (exclusive).
     */
    void write_mapping(const std::string& filename) const;

    /**
     * Replace the existing node mapping with the one loaded from the file.
     * This should be used before calling unfold(). The identifiers for new
     * duplicated nodes will follow the ones in the loaded mapping.
     */
    void read_mapping(const std::string& filename);

    /**
     * Get the id of the corresponding node in the original graph.
     */
    vg::id_t get_mapping(vg::id_t node) const;

    /**
     * Create an edge between two node orientations.
     */
    static edge_t make_edge(const HandleGraph& graph, gbwt::node_type from, gbwt::node_type to) {
        return make_pair(graph.get_handle(gbwt::Node::id(from), gbwt::Node::is_reverse(from)),
                         graph.get_handle(gbwt::Node::id(to), gbwt::Node::is_reverse(to)));
    }

private:
    /**
     * Generate a complement graph consisting of the edges that are in the
     * GBWT index but not in the input graph. Split the complement into
     * disjoint components and return the components.
     */
    std::list<bdsg::HashGraph> complement_components(MutableHandleGraph& graph, bool show_progress);

    /**
     * Generate all border-to-border paths in the component supported by the
     * indexes. Unfold the paths by duplicating the inner nodes so that the
     * paths become disjoint, except for their shared prefixes/suffixes.
     */
    size_t unfold_component(MutableHandleGraph& component, MutableHandleGraph& graph, MutableHandleGraph& unfolded);

    /**
     * Generate all paths supported by the XG index passing through the given
     * node until the border or until the path ends. Insert the generated
     * paths into the set in the canonical orientation, and use them as
     * reference paths for extending threads.
     */
    void generate_paths(MutableHandleGraph& component, vg::id_t from);

   /**
    * Generate all paths supported by the GBWT index from the given node until
    * the border. Extend paths that start/end at internal nodes using the
    * reference paths. If the node is a border node, consider all threads
    * passing through it. Otherwise consider only the threads starting from
    * it, and do not output threads reaching a border.
    */
    void generate_threads(MutableHandleGraph& component, vg::id_t from);

    /**
     * Create or extend the state with the given node orientation, and insert
     * it into the stack if it is supported by the GBWT index. Use 'starting'
     * to determine whether the initial state is for the threads starting at
     * the node or for the threads passing through the node.
     */
    void create_state(vg::id_t node, bool is_reverse, bool starting);
    bool extend_state(state_type state, vg::id_t node, bool is_reverse);

    /**
     * Try to extend the path at both ends until the border by using the
     * reference paths. Insert the extended path into the set in the canonical
     * orientation.
     */
    void extend_path(const path_type& path);

    /// Insert the path into the set in the canonical orientation.
    void insert_path(const path_type& path, bool from_border, bool to_border);

    /// Get the id for the duplicate of 'node' after 'from'.
    gbwt::node_type get_prefix(gbwt::node_type from, gbwt::node_type node);

    /// Get the id for the duplicate of 'node' before 'to'.
    gbwt::node_type get_suffix(gbwt::node_type node, gbwt::node_type to);

    /// XG and GBWT indexes for the original graph.
    const PathHandleGraph& path_graph;
    const gbwt::GBWT& gbwt_index;

    /// Mapping from duplicated nodes to original ids.
    gcsa::NodeMapping mapping;

    /// Internal data structures for the current component.
    hash_set<vg::id_t>     border;
    std::stack<state_type> states;
    std::vector<path_type> reference_paths;

    /// Tries for the unfolded prefixes and reverse suffixes.
    /// prefixes[(from, to)] is the mapping for to, and
    /// suffixes[(from, to)] is the mapping for from.
    pair_hash_map<std::pair<gbwt::node_type, gbwt::node_type>, gbwt::node_type> prefixes, suffixes;
    pair_hash_set<std::pair<gbwt::node_type, gbwt::node_type>> crossing_edges;
};

}

#endif 
