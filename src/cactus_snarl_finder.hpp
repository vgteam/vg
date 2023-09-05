///
///  \file cactus_snarl_finder.hpp
///
///  Defines a widget for finding snarls using the pinchesAndCacti library.
///

#ifndef VG_CACTUS_SNARL_FINDER_HPP_INCLUDED
#define VG_CACTUS_SNARL_FINDER_HPP_INCLUDED

#include "snarls.hpp"

namespace vg {

using namespace std;


/**
 * Class for finding all snarls using the base-level Cactus snarl decomposition
 * interface.
 */
class CactusSnarlFinder : public SnarlFinder {

protected:
    /// Holds the vg graph we are looking for sites in.
    const PathHandleGraph* graph;
    
    /// Holds the names of reference path hints
    unordered_set<string> hint_paths;
    
    /// Create a snarl in the given SnarlManager with the given start and end,
    /// containing the given child snarls in the list of chains of children and
    /// the given list of unary children. Recursively creates snarls in the
    /// SnarlManager for the children. Returns a pointer to the finished snarl
    /// in the SnarlManager. Start and end may be empty visits, in which case no
    /// snarl is created, all the child chains are added as root chains, and
    /// null is returned. If parent_start and parent_end are empty Visits, no
    /// parent() is added to the produced snarl.
    const Snarl* recursively_emit_snarls(const Visit& start, const Visit& end,
        const Visit& parent_start, const Visit& parent_end,
        stList* chains_list, stList* unary_snarls_list, SnarlManager& destination);

    /**
     * Find all the snarls with Cactus, and put them into a SnarlManager.
     * Skip breaking into connected components if "known_single_component" is true
     * Skip making the snarl manager index if finish_index is false
     */
    virtual SnarlManager find_snarls_impl(bool known_single_component, bool finish_index);
    
public:
    /**
     * Make a new CactusSnarlFinder to find snarls in the given graph.
     * We can't filter trivial bubbles because that would break our chains.
     *
     * Optionally takes a hint path name.
     */
    CactusSnarlFinder(const PathHandleGraph& graph, const string& hint_path = "");
        
    /**
     * Find all the snarls with Cactus, and put them into a SnarlManager.
     */
    virtual SnarlManager find_snarls();

    /**
     * Find all the snarls of weakly connected components in parallel.
     * Even single-threaded, this may be worth using as it will use less
     * memory by only considering each component in the context of itself.
     */
    virtual SnarlManager find_snarls_parallel();
    
};

}

#endif
