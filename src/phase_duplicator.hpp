#ifndef VG_PHASE_DUPLICATOR_H
#define VG_PHASE_DUPLICATOR_H

/** \file
 *
 * This file contains the PhaseDuplicator, which duplicates and un-crosslinks
 * complex regions of a graph using information about phased traversals.
 */

#include <utility> 
#include <vector>

#include <xg.hpp>

#include "vg.pb.h"
#include "types.hpp"

namespace vg {

using namespace std;

/**
 * Transforms complex subregions of a graph into collections of disconnected
 * distinct traversal haplotypes. Requires an xg index with a gPBWT in order to
 * work.
 */
class PhaseDuplicator {
public:
    /// Make a new PhaseDuplicator backed by the gievn index with gPBWT
    PhaseDuplicator(const xg::XG& index);
    
    /**
     * Duplicate out the subgraph induced by the given set of node IDs. Border
     * nodes of the subgraph, which edges out to the rest of the graph will be
     * identified, and all unique traversals from one border node to another
     * will be generated as new nodes, with edges connecting them to material
     * outside the graph and Translations embedding them in the original graph.
     *
     * TODO: also generate one-end-anchored traversals and internal not-
     * attached-to-a-border traversals, because there may be phase breaks in the
     * region.
     *
     * New IDs will be generated startign with next_id, and next_id will be
     * updated to the next ID after the IDs of all generated material.
     * 
     * TODO: invent some kind of thread safe ID allocator so we can do multiple
     * subgraphs in parallel without renumbering later.
     */
    pair<Graph, vector<Translation>> duplicate(const set<id_t>& subgraph, id_t& next_id) const;
    
    /**
     * List all the distinct haplotypes within a subgraph and their counts.
     * Reports each haplotype in only one direction.
     */
    vector<pair<xg::XG::thread_t, int>> list_haplotypes(const set<id_t>& subgraph) const;
    
    /**
     * List all the distinct haplotypes going through the given node in the given
     * orientation, and staying within the given subgraph, along with their
     * counts. Some may be suffixes of others.
     */
    vector<pair<xg::XG::thread_t, int>> list_haplotypes_through(xg::XG::ThreadMapping start_node, const set<id_t>& subgraph) const;
    
    /**
     * List all the distinct haplotypes actually beginning at the given node in
     * the given orientation, and staying within the given subgraph, along with
     * their counts. Haplotypes that begin and end at the same side may or may
     * not be reported multiple times, in differing orientations.
     */
    vector<pair<xg::XG::thread_t, int>> list_haplotypes_from(xg::XG::ThreadMapping start_node, const set<id_t>& subgraph) const;
    
    /**
     * List all the distinct haplotypes beginning at the given node, using the
     * given starting search state, that traverse through the subgraph.
     */
    vector<pair<xg::XG::thread_t, int>> list_haplotypes(xg::XG::ThreadMapping start_node,
        xg::XG::ThreadSearchState start_state, const set<id_t>& subgraph) const;
        
    /**
     * Get the traversals that represent the borders of the subgraph, through
     * which we can enter the subgraph.
     */
    set<xg::XG::ThreadMapping> find_borders(const set<id_t> subgraph) const;
    
    /**
     * Find the edges on the given side of the given oriented ThreadMapping that
     * cross the border of the given subgraph.
     */
    vector<Edge> find_border_edges(xg::XG::ThreadMapping mapping, bool on_start, const set<id_t>& subgraph) const;
    
    /**
     * Follow the given edge from the given node visited in the given
     * orientation, and return the node and orientation we land in.
     *
     * TODO: I think I may have written this logic already with one of the other
     * two oriented node types (NodeTraversal and NodeSide).
     */
    static xg::XG::ThreadMapping traverse_edge(const Edge& e, const xg::XG::ThreadMapping& prev);
    
    /**
     * Return a copy of the given thread in a canonical orientation (either
     * forward or reverse, whichever compares smaller).
     */
    static xg::XG::thread_t canonicalize(const xg::XG::thread_t& thread);
    
private:
    /// What XG index describes the graph we operate on?
    const xg::XG& index;
};

}

#endif 
