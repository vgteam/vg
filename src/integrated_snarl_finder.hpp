///
///  \file integrated_snarl_finder.hpp
///
///  Defines a widget for finding snarls using an internal implementation of Cactus graphs over HandleGraphs
///

#ifndef VG_INTEGRATED_SNARL_FINDER_HPP_INCLUDED
#define VG_INTEGRATED_SNARL_FINDER_HPP_INCLUDED

#include "snarls.hpp"

#include <functional>

namespace vg {

using namespace std;


/**
 * Class for finding all snarls using an integrated Cactus graph construction
 * algorithm.
 *
 * Does not produce any unary snarls. May leave edges in the root snarl, at the
 * ends of top-level chains.
 *
 * Does not (yet) use paths for rooting. Roots the decomposition at the simple
 * cycle or bridge tree path with the most bases of fixed sequence.
 */
class IntegratedSnarlFinder : public SnarlFinder {
protected:
    /// Holds the base graph we are looking for sites in.
    const PathHandleGraph* graph;
    
private:
    // Forward-declare this member type we use inside some functions.
    
    /// Represents a graph that starts as the graph of adjacency components in
    /// a HandleGraph, and which can be further merged. Can be used to
    /// represent a cactus graph or a bridge forest.
    ///
    /// Represents a graph of "components". Each component contains some number
    /// of handles from the backing graph, all reading into the component. Each
    /// handle connects that component to another component: the component that
    /// contains the flipped version of the handle. Each component is
    /// identified by a "head" handle.
    class MergedAdjacencyGraph;
    
public:
    /**
     * Make a new IntegratedSnarlFinder to find snarls in the given graph.
     */
    IntegratedSnarlFinder(const PathHandleGraph& graph);
    
    /**
     * Visit all snarls, including trivial snarls.
     *
     * Visits children before their parents.
     *
     * Calls the iteratee with the parent's boundaries if any, and the snarl's boundaries.
     *
     * Start handles are inward facing and end handles are outward facing.
     */
    void for_each_snarl_including_trivial_postorder_with_parent(const function<void(const pair<handle_t, handle_t>*, const pair<handle_t, handle_t>&)>& iteratee) const;
        
    /**
     * Find all the snarls, and put them into a SnarlManager. Make sure to
     * include trivial snarls to keep chains intact.
     */
    virtual SnarlManager find_snarls();
    
};

}

#endif
