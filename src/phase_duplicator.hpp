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
    pair<Graph, vector<Translation>>&& duplicate(set<id_t> subgraph, id_t& next_id) const;
    
private:
    /// What XG index describes the graph we operate on?
    const xg::XG& index;
};

}

#endif 
