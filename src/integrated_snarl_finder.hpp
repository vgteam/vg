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
class IntegratedSnarlFinder : public HandleGraphSnarlFinder {
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
     * Visit all snarls and chains, including trivial snarls and single-node
     * trivial chains.
     *
     * Calls begin_chain and end_chain when entrering and exiting chains in the
     * traversal. Within each chain, calls begin_snarl and end_snarl when
     * entering and exiting each snarl, in order. The caller is intended to
     * maintain its own stack to match up begin and end events.
     *
     * Each begin/end call receives the handle reading into/out of the snarl or
     * chain. 
     *
     * Both trivial and cyclic chains have the in and out handles the same.
     * They are distinguished by context; trivial chains have no shild snarls,
     * while cyclic chains do.
     *
     * Roots the decomposition at a global snarl with no bounding nodes, for
     * which begin_snarl is not called. So the first call will be begin_chain.
     *
     * Start handles are inward facing and end handles are outward facing.
     * Snarls must be oriented forward in their chains.
     */
    void traverse_decomposition(const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
        const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl) const;
};

}

#endif
